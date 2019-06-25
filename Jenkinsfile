pipeline {
    agent { 
        docker {
            image "angusgibson/fluidity:disco-a03b"
            label 'dockerhost'
        } 
    }
    environment {
        MPLBACKEND = 'PS'
        OMPI_MCA_btl = '^openib'
    }
    stages {
        stage('Configuring') {   
            steps { 
                sh 'FCFLAGS=-I/usr/lib/x86_64-linux-gnu/fortran/gfortran-mod-15/zoltan ./configure --enable-2d-adaptivity' 
            }
        }    
        stage('Building') {       
            steps { 
                sh 'make -j' ;
                sh 'make -j fltools' ;
                sh 'make manual'
            }
        }
        stage('Testing') {       
            steps { 
                sh 'make unittest' ;
                sh 'make THREADS=8 test' ;
                sh 'make THREADS=8 mediumtest'
                junit 'tests/test_result*xml'
            }
        }
    }
}
