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
        GFORTRAN_UNBUFFERED_ALL='1'
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
                sh 'make test' ;
                sh 'make mediumtest'
            }
        }
    }
    post {
         always {
            try {
               junit 'tests/test_result*xml'
           } catch(err) {
               step([$class: 'JUnitResultArchiver', testResults: 'tests/test_result*.xml'])
               if (currentBuild.result == 'UNSTABLE')
                   currentBuild.result = 'FAILURE'
               throw err
          }
             archiveArtifacts '*test_failures.tar.gz'
       }
    }
}
