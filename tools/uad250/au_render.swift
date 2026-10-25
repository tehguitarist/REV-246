// Headless offline renderer for the UAD EMT 250 AU (out-of-process: in-process instantiation segfaults).
//
//   au_render IN.wav OUT.wav [--rate 48000] [--tail 10] [--block 512]
//             [--set Param=Value ...] [--auto events.json] [--mono-in]
//
// Params are set by name (see au_probe) before the first block. Value is a number, or for list
// parameters the display string (e.g. Time="2.5 s"). --auto takes a JSON list of
// {"t": seconds, "param": name, "value": v} or {"t": seconds, "action": "reset"}; blocks are split at event
// times so a change lands on its sample.
import AVFoundation
import AudioToolbox

func fcc(_ s: String) -> FourCharCode { var v: UInt32 = 0; for c in s.utf8 { v = (v << 8) | UInt32(c) }; return v }
func die(_ m: String) -> Never { FileHandle.standardError.write((m + "\n").data(using: .utf8)!); exit(1) }

var args = Array(CommandLine.arguments.dropFirst())
guard args.count >= 2 else { die("usage: au_render IN.wav OUT.wav [--rate R] [--tail S] [--block N] [--set P=V] [--auto f.json] [--mono-in]") }
let inPath = args.removeFirst(), outPath = args.removeFirst()
var rate = 48000.0, tail = 10.0, block = 512, autoPath: String? = nil, monoIn = false, offline = true, realtime = false, oop = false
var sets: [(String, String)] = []
while !args.isEmpty {
    let a = args.removeFirst()
    switch a {
    case "--rate": rate = Double(args.removeFirst())!
    case "--tail": tail = Double(args.removeFirst())!
    case "--block": block = Int(args.removeFirst())!
    case "--auto": autoPath = args.removeFirst()
    case "--mono-in": monoIn = true
    case "--oop": oop = true  // out-of-process (sandboxed host: cannot reach the UAD DSP, outputs noise only)
    case "--no-offline": offline = false
    case "--realtime": realtime = true  // pace to the wall clock and stamp host time (DSP plugins may need it)
    case "--set":
        let kv = args.removeFirst().split(separator: "=", maxSplits: 1).map(String.init)
        guard kv.count == 2 else { die("--set needs Param=Value") }
        sets.append((kv[0], kv[1]))
    default: die("unknown arg \(a)")
    }
}

// Input: resample-free, must already be at --rate. Mono input is fed to both channels.
let inFile = try AVAudioFile(forReading: URL(fileURLWithPath: inPath))
guard inFile.processingFormat.sampleRate == rate else { die("input is \(inFile.processingFormat.sampleRate) Hz, --rate \(rate)") }
let nIn = Int(inFile.length)
let inBuf = AVAudioPCMBuffer(pcmFormat: inFile.processingFormat, frameCapacity: AVAudioFrameCount(nIn))!
try inFile.read(into: inBuf)
let inCh = Int(inFile.processingFormat.channelCount)
func inSample(_ ch: Int, _ i: Int) -> Float {
    if i >= nIn { return 0 }
    let c = (monoIn || inCh == 1) ? 0 : min(ch, inCh - 1)
    return inBuf.floatChannelData![c][i]
}
let nOut = nIn + Int(tail * rate)

var desc = AudioComponentDescription(componentType: kAudioUnitType_Effect, componentSubType: fcc("38AU"),
                                     componentManufacturer: fcc("!UAD"), componentFlags: 0, componentFlagsMask: 0)
let sem = DispatchSemaphore(value: 0)
var unit: AVAudioUnit?
AVAudioUnit.instantiate(with: desc, options: oop ? [.loadOutOfProcess] : []) { u, e in
    if let e = e { FileHandle.standardError.write("instantiate: \(e)\n".data(using: .utf8)!) }
    unit = u; sem.signal()
}
sem.wait()
guard let avau = unit else { die("could not instantiate UAD EMT 250") }
let au = avau.auAudioUnit
let fmt = AVAudioFormat(commonFormat: .pcmFormatFloat32, sampleRate: rate, channels: 2, interleaved: false)!
var pos = 0
// Input comes from a source node that reads the file at the engine's own position.
let src = AVAudioSourceNode(format: fmt) { _, _, frames, abl in
    let l = UnsafeMutableAudioBufferListPointer(abl)
    for ch in 0..<min(2, l.count) {
        let p = l[ch].mData!.assumingMemoryBound(to: Float.self)
        for i in 0..<Int(frames) { p[i] = inSample(ch, pos + i) }
    }
    return noErr
}
let engine = AVAudioEngine()
engine.attach(avau); engine.attach(src)
engine.connect(src, to: avau, format: fmt)
engine.connect(avau, to: engine.mainMixerNode, format: fmt)
try engine.enableManualRenderingMode(.offline, format: fmt, maximumFrameCount: AVAudioFrameCount(block))
au.isRenderingOffline = offline
try engine.start()
print("reported latency: \(au.latency * rate) samples")

guard let tree = au.parameterTree else { die("no parameter tree") }
func param(_ name: String) -> AUParameter {
    guard let p = tree.allParameters.first(where: { $0.displayName == name || $0.identifier == name }) else { die("no parameter \(name)") }
    return p
}
func setParam(_ name: String, _ v: Any) {
    let p = param(name)
    if let d = v as? Double { p.value = AUValue(d); return }
    let s = "\(v)"
    if let d = Double(s) { p.value = AUValue(d); return }
    if let i = p.valueStrings?.firstIndex(where: { $0.trimmingCharacters(in: .whitespaces) == s }) { p.value = AUValue(i); return }
    die("bad value \(s) for \(name)")
}
for (k, v) in sets { setParam(k, v) }

struct Ev { var s: Int; var param: String?; var value: Any?; var action: String? }
var events: [Ev] = []
if let ap = autoPath {
    let j = try JSONSerialization.jsonObject(with: Data(contentsOf: URL(fileURLWithPath: ap))) as! [[String: Any]]
    for e in j {
        events.append(Ev(s: Int((e["t"] as! Double) * rate), param: e["param"] as? String, value: e["value"], action: e["action"] as? String))
    }
    events.sort { $0.s < $1.s }
}

var outFile: AVAudioFile? = try AVAudioFile(forWriting: URL(fileURLWithPath: outPath),
    settings: [AVFormatIDKey: kAudioFormatLinearPCM, AVSampleRateKey: rate, AVNumberOfChannelsKey: 2,
               AVLinearPCMBitDepthKey: 32, AVLinearPCMIsFloatKey: true, AVLinearPCMIsNonInterleaved: false],
    commonFormat: .pcmFormatFloat32, interleaved: false)
let outBuf = AVAudioPCMBuffer(pcmFormat: fmt, frameCapacity: AVAudioFrameCount(block))!

var evi = 0
var tb = mach_timebase_info_data_t(); mach_timebase_info(&tb)
let hostHz = UInt64(1e9) * UInt64(tb.denom) / UInt64(tb.numer)
let t0 = mach_absolute_time()
while pos < nOut {
    while evi < events.count && events[evi].s <= pos {
        let e = events[evi]; evi += 1
        if e.action == "reset" { au.reset() } else if let p = e.param, let v = e.value { setParam(p, v) }
    }
    var n = min(block, nOut - pos)
    if evi < events.count { n = min(n, max(1, events[evi].s - pos)) }
    if realtime {
        let target = t0 + UInt64(Double(pos + n) / rate * Double(hostHz))
        while mach_absolute_time() < target { usleep(200) }
    }
    let st = try engine.renderOffline(AVAudioFrameCount(n), to: outBuf)
    if st != .success { die("render status \(st) at sample \(pos)") }
    try outFile!.write(from: outBuf)
    pos += n
}
outFile = nil  // closes the file and finalises the header
engine.stop()
print("rendered \(nOut) samples @ \(rate) Hz -> \(outPath)")
fflush(stdout)
_exit(0)  // the plugin segfaults in its atexit handlers
