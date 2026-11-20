# Computational protein design skills

Portable **Agent Skills** (plain markdown) for proposing, critiquing, and running realistic computational protein design campaigns.

This is **not** a Cursor-only pack. Each skill is a folder with a `SKILL.md` that any LLM agent can load: Claude Code, Codex, Cursor, Copilot, or a custom agent.

```
skills/
  protein-design-campaign/   # end-to-end campaign plans
  rfd4-proteina/             # RFD4-Proteina generation, YAML, search
  protean-editing/           # text-guided 3D edits (PROTEAN)
  momst-multiobjective/      # Pareto / test-time scaling
  protein-design-filters/    # wet-lab gates
  inverse-folding-classics/  # ProteinMPNN, AF, ESMFold, Boltz
  structure-free-design/     # Talus / Ptarmigan-1
```

Format: [Agent Skills](https://agentskills.io) — YAML frontmatter (`name`, `description`) plus instructions. Extra files (`yaml-recipes.md`, `campaign-template.md`) are progressive disclosure.

## Use in a public GitHub repo

1. Copy the `skills/` directory (and optionally `AGENTS.md`) into your repository.
2. Point agents at it. `AGENTS.md` in this folder already says: read `skills/*/SKILL.md` for protein design work.
3. Commit and push. No Cursor config is required.

Vendoring example:

```bash
git clone https://github.com/<you>/<this-repo>.git
cp -R <this-repo>/skills your-repo/skills
cp <this-repo>/AGENTS.md your-repo/AGENTS.md
```

## Load in common agents

| Agent | How |
|---|---|
| Any / GitHub | Keep `skills/` in the repo; mention it in `AGENTS.md` or the PR/issue prompt |
| Claude Code | Copy or symlink `skills/*` → `.claude/skills/` |
| Codex | Copy or symlink `skills/*` → `.agents/skills/` |
| Cursor | Copy or symlink `skills/*` → `.cursor/skills/` |
| Skills CLI | `npx skills add <org>/<repo> --skill protein-design-campaign` (after you publish) |

Symlink from a clone (one location, many adapters):

```bash
# from your repo root
ln -s ../../skills/protein-design-campaign .claude/skills/protein-design-campaign
```

## Skill map

| User asks to… | Load |
|---|---|
| Plan or critique a full campaign | `skills/protein-design-campaign/SKILL.md` |
| Run RFD4-Proteina / YAML / partial diffusion / search | `skills/rfd4-proteina/SKILL.md` |
| Local 3D edit from natural language | `skills/protean-editing/SKILL.md` |
| Balance stability / solubility / function | `skills/momst-multiobjective/SKILL.md` |
| Rank designs before ordering DNA | `skills/protein-design-filters/SKILL.md` |
| Inverse fold / classic RosettaCommons loop | `skills/inverse-folding-classics/SKILL.md` |
| IDPs, cryptic sites, ligand screens without a fold | `skills/structure-free-design/SKILL.md` |

Sources: [RFD4-Proteina](https://rfd4-proteina.org/), [PROTEAN](https://arxiv.org/abs/2610.17765), [MoMST](https://github.com/MingYangi/MoMST), [Talus Bio](https://www.talus.bio/).
