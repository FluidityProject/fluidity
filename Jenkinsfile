pipeline {
    agent { 
        docker {
            image "fluidity/baseimages:bionic-python3"
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
                sh './configure --enable-2d-adaptivity' 
            }
        }    
        stage('Building') {       
            steps { 
                sh 'make -j' ;
                sh 'make -j fltools' ;
                sh 'make makefiles' ;
                sh 'test -z "$(git status --porcelain */Makefile.dependencies)"' ;
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
