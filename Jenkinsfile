def PROJECT_NAME=JOB_NAME.split("\\s|/")[2].toLowerCase()

pipeline {
    agent { 
        docker {
            image "fluidity/baseimages:${PROJECT_NAME}"
            label 'azure-linux'
        } 
    }
    environment {
        MPLBACKEND = 'PS'
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
                sh 'make -j 8 unittest' ;
                sh 'make -j 8 test' ;
                sh 'make -j 8 mediumtest'
            }
        }
    }
    post {
        always {
            junit 'tests/test_result*xml'
        }
    }
}
