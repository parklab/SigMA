file_bed_list=$1

for file in `cat $file_bed_list`
  do 
    dir=`echo "$(dirname "${file}")"`
    cd $dir
    file_name=`echo "$(basename "${file}")"`
    git lfs pull --include=$file_name
    cd -
  done
  