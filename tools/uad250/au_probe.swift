// Headless AU probe: instantiate UAD EMT 250, list parameters, run a short render.
import AVFoundation
import AudioToolbox

var desc = AudioComponentDescription(componentType: kAudioUnitType_Effect,
    componentSubType: 0x38415520 >> 0 == 0 ? 0 : FourCharCode("38AU"),
    componentManufacturer: FourCharCode("!UAD"), componentFlags: 0, componentFlagsMask: 0)

extension FourCharCode {
    init(_ s: String) { var v: UInt32 = 0; for c in s.utf8 { v = (v << 8) | UInt32(c) }; self = v }
}
desc.componentSubType = FourCharCode("38AU")

let sem = DispatchSemaphore(value: 0)
var unit: AVAudioUnit?
AVAudioUnit.instantiate(with: desc, options: CommandLine.arguments.contains("oop") ? [.loadOutOfProcess] : []) { u, e in
    unit = u
    if let e = e { print("instantiate error:", e) }
    sem.signal()
}
sem.wait()
guard let au = unit else { print("no unit"); exit(1) }
print("name:", au.name, "latency(s):", au.auAudioUnit.latency, "tail(s):", au.auAudioUnit.tailTime)
print("maxFrames:", au.auAudioUnit.maximumFramesToRender)
if let tree = au.auAudioUnit.parameterTree {
    for p in tree.allParameters {
        print(String(format: "%llu\t%@\t[%g..%g]\tval=%g\tunit=%d\tvals=%@", p.address, p.displayName, p.minValue, p.maxValue, p.value, p.unit.rawValue, (p.valueStrings ?? []).joined(separator: "|")))
    }
}
print("presets:", au.auAudioUnit.factoryPresets?.map { "\($0.number):\($0.name)" } ?? [])
