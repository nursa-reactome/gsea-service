Reactome GSEA Service
=====================
The Reactome GSEA service performs
[GSEA](http://software.broadinstitute.org/gsea/index.jsp)
enrichment analysis on an expression dataset against the
[Reactome Pathways Gene Set](https://reactome.org/download/current/ReactomePathways.gmt.zip).

This is a REST service with the _analyse_` `POST` method. The request body
is a list consisting of [_symbol_, _value_] records, where:

* _symbol_ is the gene symbol

* _value_ is the rank value, e.g. fold change or p-value

Optional query parameters include the following:

* _nperms_: number of permutations

* _dataSetSizeMin_: minimum gene set size

* _dataSetSizeMax_: maximum gene set size

The procedure and parameters are described in the GSEA User Guide
[Running Analyses](http://software.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Running_Analyses)
topic.

Installation
------------
1. Clone this Git repository.

2. Install [Maven](http://maven.apache.org/install.html).

3. Build the `.war` file:

        mvn clean package -U

6. Copy the `.war` file to tomcat:

        cp target/GseaService.war $TOMCAT_HOME/webapps

   where `$TOMCAT_HOME` is the tomcat deployment location.

7. Alternatively, the REST service can be started locally in an embedded
   server with the Maven `spring-boot:run` goal:

        mvn tomcat7:run
