#!/usr/bin/env python3
"""
rss_monitor.py — Background RSS peak tracker for build_2step_pipeline.sh

Polls /proc every INTERVAL seconds, sums RSS across the entire process
tree rooted at WATCH_PID (the shell script), and records the peak.

The shell script writes the current step name to <peak_file>.label so
the monitor can report WHICH step hit the peak.

Usage (called by build_2step_pipeline.sh — not meant to be run directly):
    python3 rss_monitor.py <watch_pid> <peak_file> [interval_sec]
"""
import sys, os, time

def get_all_pids():
    pids = []
    try:
        for e in os.listdir('/proc'):
            if e.isdigit(): pids.append(int(e))
    except Exception: pass
    return pids

def get_ppid(pid):
    try:
        with open(f'/proc/{pid}/status') as f:
            for line in f:
                if line.startswith('PPid:'):
                    return int(line.split()[1])
    except Exception: pass
    return -1

def get_rss_kb(pid):
    try:
        with open(f'/proc/{pid}/status') as f:
            for line in f:
                if line.startswith('VmRSS:'):
                    return int(line.split()[1])
    except Exception: pass
    return 0

def get_tree_rss(root_pid):
    all_pids = get_all_pids()
    parent = {}
    for p in all_pids:
        pp = get_ppid(p)
        if pp >= 0:
            parent[p] = pp
    tree = set()
    def collect(pid):
        tree.add(pid)
        for p, pp in parent.items():
            if pp == pid and p not in tree:
                collect(p)
    collect(root_pid)
    return sum(get_rss_kb(p) for p in tree)

def main():
    if len(sys.argv) < 3:
        sys.exit(1)
    watch_pid = int(sys.argv[1])
    peak_file = sys.argv[2]
    interval  = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
    label_file = peak_file + ".label"
    peak_kb = 0
    peak_label = "unknown"

    while os.path.exists(f'/proc/{watch_pid}'):
        rss_kb = get_tree_rss(watch_pid)
        try:
            with open(label_file) as f: lbl = f.read().strip()
        except Exception: lbl = peak_label

        if rss_kb > peak_kb:
            peak_kb   = rss_kb
            peak_label = lbl
            try:
                with open(peak_file, 'w') as f:
                    f.write(f'{peak_kb}\t{peak_label}\n')
            except Exception: pass

        time.sleep(interval)

if __name__ == '__main__':
    main()
