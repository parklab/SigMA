usage="$(basename "$0") [bed_list] [mut_file] [output_file]"

while getopts ':hbio:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
  esac
done

#shift $((OPTIND - 1))
bed_file_list=$1
input_file=$2
output_file=$3

echo $bed_file_list
echo $input_file
echo $output_file

if [ -f $output_file ]; then
  echo "Output file $file already exists. Use a different name."
  exit
fi

cat $bed_file_list | while read bed_file
do 
  echo $bed_file
  bedtools intersect -a $input_file -b $bed_file >> $output_file
done

Rscript make_table_indel_repeat.R $output_file counts_$output_file.csv
