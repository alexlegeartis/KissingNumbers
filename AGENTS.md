# Agent guidelines

Leave at least one logical CPU free from agent compute workloads. On the current
128-CPU host, reserve CPU 127 and pin compute commands to CPUs 0–126:

```sh
taskset -c 0-126 <command>
```

Apply this affinity to every concurrent compute job, including GPU jobs with CPU
workers; child processes inherit it. Keep the combined worker and thread budget
at or below 127, accounting for nested BLAS/OpenMP thread pools. On another host,
check the available CPU affinity and reserve at least one available logical CPU.

Use the NVIDIA RTX A4000 with 16 GB VRAM for GPU compute. CUDA access to this GPU
was verified on this host on 2026-11-04 outside the execution sandbox.
Inside the sandbox, GPU device access was unavailable. Use the execution tool's
supported escalated mode when GPU access requires it, within the configured
permissions. Recheck GPU availability and device indices before launching work
on a new host or session.
