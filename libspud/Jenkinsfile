pipeline {

    agent { 
        docker {
	    image 'fluidity/spudbase:bionic'
            label 'dockerhost'
        } 
    }
    environment {
        OMPI_MCA_btl = '^openib'
    }
    stages {
        stage('Configuring') {   
            steps { 
                sh './configure --prefix=/home/spud' 
            }
        }    
        stage('Building (python3 default)') {       
            steps { 
                sh 'make -j'
		sh 'make install'
            }
        }
	stage('Testing fortran library') {
            steps { 
                sh 'make junittest' ;
                junit 'src/tests/test_result*xml'
            }
	}
	stage('Testing (python3 default)') {       
            steps { 
                withEnv(['PYTHONPATH=/home/spud/lib/python3.6/site-packages',
		         'LD_LIBRARY_PATH=/home/spud/lib']) {
		    sh 'cd python; python3 test_libspud_junit.py'
                }
                junit 'python/test_result*xml'
		sh 'rm python/test_result*xml'
            }
        }
	stage('Building libspud for python2') {       
            steps { 
		sh 'cd python; python2 setup.py install --prefix=/home/spud; cd ..'
            }
        }
        stage('Testing for python2') {       
            steps { 
	        withEnv(['PYTHONPATH=/home/spud/lib/python2.7/site-packages',
		         'LD_LIBRARY_PATH=/home/spud/lib']) {
		    sh 'cd python; python2 test_libspud_junit.py' 
                }
                junit 'python/test_result*xml'
		sh 'rm python/test_result*xml'
            }
        }
	stage('Building documentation') {       
            steps { 
		sh 'make doc'
            }
        }
    }
}
