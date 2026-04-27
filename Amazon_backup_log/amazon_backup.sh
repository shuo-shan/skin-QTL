# 1. compile table of files to backup: 
##   note: make sure to utilize the amazon_backups.xlsx file in our Dropbox folder to make this list
##   field1: file path
##   field2: file name
##   field3: s3bucket path

# 2. paste the table of files into a file here (main folder) "s3_backup_files"

# 3. run this command at the main folder
mainpth=$PWD
cd $mainpth
rm source_backup
while read line;do
  pth=$(echo $line | cut -d' ' -f1)
  f=$(echo $line | cut -d' ' -f2)
  amazonbin=$(echo $line | cut -d' ' -f3)
  i=${pth}/${f}

  cd $pth
  md5sum $f > $f.md5sum
  echo $f "done md5sum"

  cd $mainpth
  echo ${i} 
  echo "mysql -h galaxy.umassmed.edu -u amazon -pamazon2017 biocore -e 'INSERT INTO amazon_backup (file_name, s3bucket, date_created, date_modified) VALUES (\"$i\", \"$amazonbin\", now(), now())'" >> source_backup
done < s3_backup_files
source source_backup

