f=$1
dir=/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq
outputName=linecount.txt

# count line
wcl=$(zcat "$f" | wc -l)

# full path to output file
output_file=${dir}/${outputName}

# acquire lock before writing
lock_file=${output_file}.lock

# open the lock file and associate with file descriptor 9
exec 9>"$lock_file" # lock file is created if it does not already exist.

# acquire exclusive lock
flock 9 # apply an advisory lock on this file

# write to the output file
echo -e "${f}\t${wcl}" >> "$output_file"

# release the lock (optional since script ends)
flock -u 9

# Close file descriptor (optional, but tidy)
exec 9>&-

