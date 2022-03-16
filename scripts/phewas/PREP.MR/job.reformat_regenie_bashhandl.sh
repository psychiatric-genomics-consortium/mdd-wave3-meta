printColumns () 
{ 
read names
while read $names; do
    for col in $*
    do
        eval "printf $col '%s ' \$$col"
    done
    echo
done
}