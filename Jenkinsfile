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
            }
        }
        stage('Testing') {       
            steps { 
                sh './bin/testharness -f netcdf_read_errors.xml'
                sh 'cat tests/netcdf_read_errors/valid.err'
                sh './bin/testharness -f unresolvable_mesh_dependency.xml'
                sh './bin/testharness -f unresolvable_diagnostic_dependency.xml'
                sh './bin/testharness -f unresolvable_diagnostic_dependency2.xml'
                sh 'cat tests/unresolvable*/fluidity.err*'
            }
        }
    }
    post {
        always {
            junit 'tests/test_result*xml'
        }
    }
}
