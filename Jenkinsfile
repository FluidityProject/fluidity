pipeline {
    agent none
    environment {
        MPLBACKEND = 'PS'
        OMPI_MCA_btl = '^openib'
    }
    stages {
        stage('Testing') {
            parallel {
                // For each combination of parameters required, build and test
                stage('Build and test xenial autoconf container') {
                     agent { 
                        docker {
                           image "fluidity/baseimages:xenial"
                           label 'dockerhost'
                        }
                     } 
                     steps {
                        slackSend "Build started - ${env.JOB_NAME} ${env.BUILD_NUMBER} (<${env.BUILD_URL}|Open>)"
                        sh './configure --enable-2d-adaptivity' 
                        sh 'make' ;
                        sh 'make fltools' ;
                        sh 'make manual'
                        sh 'make unittest' ;
                        sh 'make THREADS=8 shorttest' ;
                        sh 'make THREADS=8 mediumtest'
                        junit 'tests/test_result*xml'
                     }
                }
                stage('Build and test xenial cmake container') {
                     agent { 
                        docker {
                           image "fluidity/baseimages:xenial"
                           label 'dockerhost'
                        }
                     } 
                     environment {
                        PETSc_DIR = '/usr/lib/petscdir/3.6.3/linux-gnu-c-opt/lib/petsc/conf'
                     }
                     steps {
                        slackSend "Build started - ${env.JOB_NAME} ${env.BUILD_NUMBER} (<${env.BUILD_URL}|Open>)"
                        sh 'cmake -D USE_LIBMBA=ON -DUSE_EXODUSII=ON .' 
                        sh 'make' ;
                        sh 'make fltools' ;
                        sh 'make unittest' ;
                        sh 'make THREADS=8 shorttest' ;
                        sh 'find . -name test_results.xml' ;
                        junit 'test_results.xml'
                     }
                }
           }
       }
    }        
    post {
        aborted {
            slackSend(color: '#DEADED',
	              message: "Build aborted - ${env.JOB_NAME} ${env.BUILD_NUMBER} (<${env.BUILD_URL}|Open>)")
        }
	success {
	    slackSend (color: 'good',
	     message: "Build completed successfully - ${env.JOB_NAME} ${env.BUILD_NUMBER} (<${env.BUILD_URL}|Open>)")
        }
	unstable {
	    slackSend(color: 'warning',
	              message: "Build completed with test failures - ${env.JOB_NAME} ${env.BUILD_NUMBER} (<${env.BUILD_URL}|Open>)")
            script {
                currentBuild.result = "FAILURE"
            }
        }
	failure {
	    slackSend(color: 'danger',
	              message: "Build failed - ${env.JOB_NAME} ${env.BUILD_NUMBER} (<${env.BUILD_URL}|Open>)")
        }
    }
}
