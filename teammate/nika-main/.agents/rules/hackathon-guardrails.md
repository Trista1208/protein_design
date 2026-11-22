---
trigger: always_on
---
1. NEVER edit, rename, or delete any existing file outside `hackathon_extensions/` without explicit user confirmation.
2. Maintain `hackathon_extensions/AGENT_PROGRESS.md` after every completed file so work can resume cleanly from saved state if interrupted.
3. Do not run long GPU jobs or long API polling loops in shell commands; only write modular code, dry-run/mock tests, and docs.