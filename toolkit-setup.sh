#!/bin/bash
if [[ $# -gt 1 ]]; then
  PARAM=$2
else 
  PARAM=$1
fi
case $PARAM in
apt ) 
    cd /tmp
    wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
    # add to your apt sources keyring so that archives signed with this key will be trusted.
    sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
    # remove the public key
    rm GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
    # add intel package repository
    sudo add-apt-repository "deb https://apt.repos.intel.com/oneapi all main"
    # install components
    while getopts ":cr" opt; do
        case $opt in
        c )
            sudo apt -y install intel-basekit intel-hpckit
            ;;
        r )
            sudo apt -y install intel-oneapi-mpi intel-oneapi-openmp
            ;;
        esac
    done
    ;;
yum ) 
    tee > /tmp/oneAPI.repo <<-EOF
    [oneAPI]
    name=Intel® oneAPI repository
    baseurl=https://yum.repos.intel.com/oneapi
    enabled=1
    gpgcheck=1
    repo_gpgcheck=1
    gpgkey=https://yum.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
EOF
    sudo mv /tmp/oneAPI.repo /etc/yum.repos.d 
    # install components
    while getopts ":cr" opt; do
        case $opt in
        c )
            sudo yum install intel-basekit intel-hpckit
            ;;
        r )
            sudo yum install intel-oneapi-mpi intel-oneapi-openmp
            ;;
        esac
    done
    ;;
dnf ) 
    tee > /tmp/oneAPI.repo <<-EOF
    [oneAPI]
    name=Intel® oneAPI repository
    baseurl=https://yum.repos.intel.com/oneapi
    enabled=1
    gpgcheck=1
    repo_gpgcheck=1
    gpgkey=https://yum.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
EOF
    sudo mv /tmp/oneAPI.repo /etc/yum.repos.d
    echo $PARAM
    # install components
    while getopts ":cr" opt; do
        case $opt in
        c )
            sudo dnf install intel-basekit intel-hpckit
            ;;
        r )
            sudo dnf install intel-oneapi-mpi intel-oneapi-openmp
            ;;
        esac
    done
    ;;
zypper )     
    sudo zypper addrepo https://yum.repos.intel.com/oneapi oneAPI
    rpm --import https://yum.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
        # install components
    while getopts ":cr" opt; do
        case $opt in
        c )
            sudo zypper install intel-basekit intel-hpckit
            ;;
        r )
            sudo zypper install intel-oneapi-mpi intel-oneapi-openmp
            ;;
        esac
    done
    ;;
* )
    echo "Unknown parameter: $PARAM"
    echo "Available parameters are: 'apt','yum','dnf','zypper'"
    ;;
esac


