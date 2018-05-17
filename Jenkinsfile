pipeline {
    agent none
    stages {
        stage('Parameterising') {
            agent { label 'azure-linux' }
            steps {
                script {
                    def version = sh script: 'git rev-parse --short HEAD | tr -d "\n"', returnStdout: true
                    def username = sh script: 'whoami | tr -d "\n"', returnStdout: true
                    def userid = sh script: 'id -u | tr -d "\n"', returnStdout: true
                    def imageTag = "${env.BUILD_NUMBER}_${version}"
                }
            }
        }
        stage('Building') {
            parallel {
                stage('Centos') {
                    agent { label 'azure-linux' }
                    steps { buildAndTest('centos') }
                }
                stage('Ubuntu') { 
                    agent { label 'azure-linux' }
                    steps { buildAndTest('ubuntu') }
                }
            }
        }
    }
}

def buildAndTest (def osBase) { 
    script {
        def imageName = "fluidity/jenkins-${osBase}"
        sh "cp docker/Dockerfile.${osBase} Dockerfile"
        def fluidityImage = docker.build("${imageName}:${imageTag}", "--build-arg userid=${userid} --build-arg username=${username} .")
        fluidityImage.inside() {
            sh './configure --enable-2d-adaptivity'
            sh 'make'
            sh 'make fltools'
            sh 'make manual'
            sh 'set +e ; make unittest'
            sh 'set +e ; ../bin/testharness -x test_results.xml -l short $(EXCLUDE_TAGS) -n $(THREADS)s '
            sh 'set +e ; make mediumtest'
        }
    }
}
