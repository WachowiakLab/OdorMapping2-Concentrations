#!/usr/bin/env bash
# Install the highconc-figure skill for Claude Code.
#
#   ./install-claude-skill.sh          -> ~/.claude/skills  (available everywhere)
#   ./install-claude-skill.sh .claude/skills   -> just this project
#
# Copies the CURRENT style files in beside the skill, so the installed copy is
# self-contained and can never be stale relative to what you installed today.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STYLE_DIR="$(dirname "$HERE")"
DEST="${1:-$HOME/.claude/skills}/highconc-figure"

mkdir -p "$DEST"
cp "$HERE/claude-skill/highconc-figure/SKILL.md" "$DEST/"
cp "$HERE/highconc_style.py"                     "$DEST/"
cp "$STYLE_DIR/STYLE_GUIDE.md"                   "$DEST/"
cp "$STYLE_DIR/CHECKLIST.md"                     "$DEST/"
mkdir -p "$DEST/palette"
cp "$STYLE_DIR"/palette/*.csv                    "$DEST/palette/"

echo "installed -> $DEST"
printf '  %s\n' SKILL.md highconc_style.py STYLE_GUIDE.md CHECKLIST.md "palette/*.csv"
echo
echo "Restart Claude Code, then ask for a figure — or type /highconc-figure"
