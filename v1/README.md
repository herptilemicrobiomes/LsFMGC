
How to take list of sample names to keep and join to the samples.csv data
Datasets need to be sorted:
I have made these file copies in the `lib` folder
```
sort wood_frog_samples.txt > s
mv s wood_frog_samples.txt
```

Then run join, it will default to using first column in each dataset as join, columns are delimited by ','

```
cd lib
join -t,  wood_frog_samples.txt  samples.csv > ../wood_frog_samples.csv
```
