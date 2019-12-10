#!/bin/bash

wget -O DOIRootCA2.cer http://sslhelp.doi.net/docs/DOIRootCA2.cer

keytool -import -file DOIRootCA2.cer -alias DOIRootCA2 -keystore cacerts -keypass changeit -storepass changeit -noprompt
printf 'org.gradle.jvmargs=-Djavax.net.ssl.trustStore=' > gradle.properties
echo -n `pwd` >> gradle.properties
printf "/cacerts -Djavax.net.ssl.trustStorePassword=changeit\n" >> gradle.properties
cat gradle.properties