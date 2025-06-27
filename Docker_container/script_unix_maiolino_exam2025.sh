if test -f "./configurationFile_maiolino_exam2025.txt"; then
echo "$FILE exists."
else
pwd > configurationFile_maiolino_exam2025.txt
fi
tt=$(head configurationFile_maiolino_exam2025.txt)
mkdir $tt
cp ./configurationFile_maiolino_exam2025.txt $tt
rm $tt/id_maiolino_exam2025.txt
docker run --platform linux/amd64 -itv $tt:/sharedFolder -v /var/run/docker.sock:/var/run/docker.sock --cidfile  $tt/id_maiolino_exam2025.txt --privileged=true -p  8888:8888 maiolinoaurelio/pab_exam2025
