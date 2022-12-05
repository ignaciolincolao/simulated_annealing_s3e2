thread=(32 256 512 1024)


for hilo in "${thread[@]}"
do
    for count in {0..9}
    do
        line=$(head -n 1 "./save/$count"_"$hilo"_"$hilo-info-test.txt")
        echo "$line" >> info_save_test.txt
    done
done