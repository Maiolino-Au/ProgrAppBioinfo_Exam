@Set "Build=%CD%"
@Echo(%Build%
@If Not Exist "configurationFile_maiolino_exam2025.txt" Set /P "=%Build%" 0<NUL 1>"configurationFile_maiolino_exam2025.txt"
mkdir %Build%
copy configurationFile_maiolino_exam2025.txt %Build%
del %Build%\id_maiolino_exam2025.txt
docker run --platform linux/amd64 -itv %Build%:/sharedFolder -v /var/run/docker.sock:/var/run/docker.sock --privileged=true --cidfile  %Build%\id_maiolino_exam2025.txt  -p  8888:8888 maiolino_exam2025
