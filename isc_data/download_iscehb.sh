# ISC-EHB data downloader
# written by Nanqiao Du

# the year of data you want to download
startyear=2016
endyear=2016


# url
url=ftp://isc-mirror.iris.washington.edu/pub/isc-ehb

# main module
for((y=$startyear;y<=$endyear;y++));do
    name=$y.res.gz
    echo $name
    wget $url/$name
    gunzip $name
done
