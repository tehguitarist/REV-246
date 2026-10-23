"""Run a command under a hard memory cap (NEXT_STEPS §2: every run that loads library data or fits anything). macOS
ignores ulimit -v and RLIMIT_AS, so this polls the physical footprint (resident + compressed + swapped, as Activity
Monitor's "Memory") of the command's whole process group, workers included, and kills the group above the cap.

    .venv/bin/python tools/memguard.py [--limit-gb 6] [--poll 0.2] -- <command ...>

Prints the peak footprint at the end (use it to measure one job before scaling up). Exit code: the command's, or 137
when the cap killed it.
"""
import argparse, ctypes, os, signal, subprocess, sys, time

_libproc = ctypes.CDLL("/usr/lib/libproc.dylib", use_errno=True)
RUSAGE_INFO_V2 = 2


class _RusageV2(ctypes.Structure):
    _fields_ = [("uuid", ctypes.c_uint8 * 16)] + [(n, ctypes.c_uint64) for n in (
        "user_time", "system_time", "pkg_idle_wkups", "interrupt_wkups", "pageins", "wired_size", "resident_size",
        "phys_footprint", "proc_start_abstime", "proc_exit_abstime", "child_user_time", "child_system_time",
        "child_pkg_idle_wkups", "child_interrupt_wkups", "child_pageins", "child_elapsed_abstime", "diskio_bytesread",
        "diskio_byteswritten")]


def footprint(pid):
    """Physical footprint of one process in bytes (0 if it has gone)."""
    r = _RusageV2()
    if _libproc.proc_pid_rusage(pid, RUSAGE_INFO_V2, ctypes.byref(r)) != 0:
        return 0
    return r.phys_footprint


def group_pids(pgid):
    out = subprocess.run(["ps", "-A", "-o", "pid=,pgid="], capture_output=True, text=True).stdout.split()
    return [int(p) for p, g in zip(out[::2], out[1::2]) if int(g) == pgid]


def main():
    ap = argparse.ArgumentParser(usage=__doc__)
    ap.add_argument("--limit-gb", type=float, default=6.0)
    ap.add_argument("--poll", type=float, default=0.2)
    ap.add_argument("cmd", nargs=argparse.REMAINDER)
    a = ap.parse_args()
    cmd = a.cmd[1:] if a.cmd[:1] == ["--"] else a.cmd
    if not cmd:
        ap.error("no command")
    limit = a.limit_gb * 2 ** 30
    p = subprocess.Popen(cmd, start_new_session=True)       # own process group: workers inherit it
    peak, killed = 0, False
    t0 = time.time()
    while p.poll() is None:
        tot = sum(footprint(q) for q in group_pids(p.pid))
        peak = max(peak, tot)
        if tot > limit:
            os.killpg(p.pid, signal.SIGKILL)
            killed = True
            print(f"[memguard] KILLED: footprint {tot / 2 ** 30:.2f} GB > cap {a.limit_gb:g} GB after {time.time() - t0:.0f} s",
                  file=sys.stderr, flush=True)
            break
        time.sleep(a.poll)
    rc = p.wait()
    try:                                                     # stragglers (workers) of a group whose leader exited
        os.killpg(p.pid, signal.SIGKILL)
    except OSError:                                          # group gone (macOS may say EPERM for an emptied group)
        pass
    print(f"[memguard] peak footprint {peak / 2 ** 30:.2f} GB, {time.time() - t0:.0f} s, exit {137 if killed else rc}",
          file=sys.stderr, flush=True)
    sys.exit(137 if killed else rc)


if __name__ == "__main__":
    main()
