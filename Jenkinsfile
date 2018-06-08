def PROJECT_NAME=JOB_NAME.split("\\s|/")[2].toLowerCase()

pipeline {
    agent { 
        docker {
            image "fluidity/baseimages:${PROJECT_NAME}"
            label 'azure-linux-8core'
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
                sh 'make manual'
            }
        }
        stage('Testing') {       
            steps { 
                sh 'make unittest' ;
                sh 'make THREADS=8 test' ;
                sh 'make THREADS=8 mediumtest'
            }
        }
    }
    post {
        always {
            junit 'tests/test_result*xml'
        }
    }
}
