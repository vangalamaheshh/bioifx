#!/bin/bash

# Moves Miseq data from Isilon to Aberdeen
# Ignores transfering .jpg files

COPY_TO_DIR="/cores/mbcf/mbcf-storage/MiSeq/";
RSYNC_DIRS="";

echo "Start time: $(date)"; >&2
echo "------------------"; >&2

for dir in "$@"; do
	dir=$( echo "$dir" | sed -e "s/\/$//" );
	RSYNC_DIRS="$RSYNC_DIRS $dir";
done

rsync --archive --progress --compress --filter='exclude *.jpg' $RSYNC_DIRS zth1@bigzfs1:$COPY_TO_DIR

echo "End time: $(date)"; >&2
echo "-------------------"; >&2

exit $?;
