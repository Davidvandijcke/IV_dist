#!/bin/bash
# Auto-sync coefficient_projection_v1.tex to Overleaf folder

SOURCE="/Users/davidvandijcke/University of Michigan Dropbox/David Van Dijcke/IV_dist/coefficient_projection_v1.tex"
DEST="/Users/davidvandijcke/University of Michigan Dropbox/David Van Dijcke/Apps/Overleaf/IV_distributions/coefficient_projection_v1.tex"

# Watch for changes and copy
/opt/homebrew/bin/fswatch -o "$SOURCE" | while read; do
    cp "$SOURCE" "$DEST"
    echo "$(date): Synced to Overleaf"
done
