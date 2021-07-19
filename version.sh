rev=$(git log --pretty=format:'%h' -n 1)
sed -i '7s/.*/    character(7), parameter :: revision = "'$rev'"/' parameters.f90
