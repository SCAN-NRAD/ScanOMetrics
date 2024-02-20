sphinx-apidoc -M --tocfile api --force -o ./source ../scanometrics
sphinx-build -b html ./source ./build
