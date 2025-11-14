coolingrates=(0.98 0.989 0.998)
temp_inital=(100000 1000 10)

mkdir "./save_example_15_30_25"
cd ./build
let seed=125
for repeat in {1..30}
do 
    for rate in "${coolingrates[@]}"
    do
        mkdir "../save_example_15_30_25/rate_$rate"
        echo "---------######## Estado 15 30 25 ####### En rate_$rate donde $repeat de 30 ######----------"
        ((seed+=$repeat*2))
# Orden de los parametros
# ejecutable | temp_ini | min_temp | alpha1 | alpha2 | alpha3 | coolingRate | k_reheating | e_const | n_reheating | len1 | len2 | len3 | len4 | Th | n_block | n_thread | rutaSave | prefijo | name_exp
        ./s3e2_sas 100000 0.00000009 15 30 25 $rate 0.0   1 1 "../save_example_15_30_25/rate_$rate/" "$repeat" $seed
    done
done
echo "---------######## Completado 15 30 25 ####### Completado rate_$rate donde $repeat de 30 ######----------"

