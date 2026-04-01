#!/usr/bin/env python3
import os
import shutil

# ----------------- SETTINGS -----------------
START_DIR   = "."                 # where to start searching (current dir)
TARGET_NAME = "frame1000.pdb"     # file to look for
DEST_DIR    = "all_frame1000"     # new parent directory to collect copies
NUM_DIRS    = 3                   # how many last subdirs to include in name
# --------------------------------------------

def make_safe(segment: str) -> str:
    """
    Make a path segment safe for filenames.
    """
    # replace spaces and slashes and other annoying chars
    bad_chars = [' ', '/', '\\', ':', ';']
    for ch in bad_chars:
        segment = segment.replace(ch, "_")
    return segment

def main():
    os.makedirs(DEST_DIR, exist_ok=True)

    copied = 0
    for root, dirs, files in os.walk(START_DIR):
        if TARGET_NAME not in files:
            continue

        src_path = os.path.join(root, TARGET_NAME)

        # Get last NUM_DIRS directory names from root
        # e.g. .../V5/23/rep1  -> ["V5", "23", "rep1"]
        parts = os.path.normpath(root).split(os.sep)
        if len(parts) >= NUM_DIRS:
            tail = parts[-NUM_DIRS:]
        else:
            tail = parts

        tail = [make_safe(p) for p in tail]
        base_name = "__".join(tail) + "_frame1000.pdb"
        dest_path = os.path.join(DEST_DIR, base_name)

        # Avoid overwriting if name collides
        if os.path.exists(dest_path):
            i = 1
            while True:
                alt_name = "__".join(tail) + f"_frame1000_{i}.pdb"
                alt_path = os.path.join(DEST_DIR, alt_name)
                if not os.path.exists(alt_path):
                    dest_path = alt_path
                    break
                i += 1

        print(f"Copying: {src_path}  ->  {dest_path}")
        shutil.copy2(src_path, dest_path)
        copied += 1

    print(f"\nDone. Copied {copied} files into '{DEST_DIR}'.")

if __name__ == "__main__":
    main()
