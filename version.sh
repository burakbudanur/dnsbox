rev=$(git log --pretty=format:'%h' -n 1)
sed -i '9s/.*/    character(7) :: revision="'$rev'"/' m_parameters.f90
