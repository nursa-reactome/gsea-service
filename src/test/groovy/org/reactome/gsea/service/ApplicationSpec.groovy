package org.reactome.gsea.service

import org.junit.Test
import org.junit.experimental.categories.Category
import org.reactome.gsea.service.IntegrationTest

import java.util.regex.Pattern
import spock.lang.Specification
import spock.lang.Shared
import groovyx.net.http.RESTClient
import groovy.json.JsonSlurper

/**
 * @author Fred Loney <loneyf@ohsu.edu>
 * 
 * The GSEA analysis test.
 * 
 * Note: this test spec is disabled by default. The pom.xml gmavenplus-plugin
 * goals were commented out by Guanming since they resulted in some error
 * in his enviroment. In order to run the tests, restore the pom.xml goals
 * but don't check in that change to the pom.xml file. 
 * 
 * Note: the pre-condition for this test is that a test gsea-server
 * is listening on port 8282.
 * Consequently, the gsea-server package must be built in Maven with
 * the skipTests flag set.
 */
// Note: the thick tangle of "convenience" Java Spring Boot testing
// annotations, e.g. @ContextConfiguration, leads to obscure errors,
// e.g. ContextConfiguration class version conflicts. The work-around
// is to avoid annotation magic and run this as a straight script
// with a test gsea-server already listening on port 8282.
@Category(IntegrationTest.class)
class ApplicationSpec extends Specification {

    final FIXTURES = "src/test/resources/fixtures"

    // The .rnk file format is <gene name><white space><rank value>
    final REGEX = /^(\w+)\s+(-?\d*(\.\d+)?(E-?\d+)?)$/

    @Shared
    def RESTClient client
    
    def setupSpec() {
        client = new RESTClient('http://localhost:8282/')
    }
    
    def "performs analysis on a small ranked list"() {
        given:
            def geneSet = "gTqItVnDEP"
            def expected = "gTqItVnDEP_no_opt"
        when:
            def response = analyse(geneSet)
        then:
            validate(geneSet, response, expected)
    }
    
    def "performs analysis using a text request"() {
        given:
            def geneSet = "gTqItVnDEP"
            def expected = "gTqItVnDEP_no_opt"
        when:
            def response = analyse_text_input(geneSet)
        then:
            validate(geneSet, response, expected)
    }

    def "applies the min dataset size option"() {
        given:
            def geneSet = "gTqItVnDEP"
            def expected = "gTqItVnDEP_min_opt"
        when:
            def response = analyse(geneSet, [dataSetSizeMin: 40])
        then:
            validate(geneSet, response, expected)
    }
    
    def analyse(basename, opts=[:]) {
        // Read the .rnk file.
        def input = this.getClass().getResource("/fixtures/${basename}.rnk")
        // Filter for the data lines.
        def data = input.readLines().findAll { it ==~ REGEX }
        // The [[name, value], ...] array.
        def records = data.collect { line ->
            def match = line =~ REGEX
            [match[0][1], match[0][2]]
        }
        // The POST JSON body.
        def json = groovy.json.JsonOutput.toJson(records)
        // Perform the GSEA analysis.
        return client.post(
            path: "analyse",
            query: opts,
            body: json,
            contentType: "application/json",
            requestContentType: "application/json"
        )
    }
    
    def analyse_text_input(basename, opts=[:]) {
        // Read the .rnk file into a string.
        def payload = new File(FIXTURES + "/${basename}.rnk").text
        // Perform the GSEA analysis.
        return client.post(
            path: "analyse",
            query: opts,
            body: payload,
            contentType: "application/json",
            requestContentType: "text/plain"
        )
    }
    
    def validate(basename, response, fixture) {
        // The expected input URL.
        def url = this.getClass().getResource("/fixtures/${fixture}.json")
        // The expected data as a set of Map objects.
        def expected = new groovy.json.JsonSlurper().parseText(url.text).collect {
            [pathway: it.pathway, regulationType: it.regulationType]
        } as Set
                    
        // Check the response status code.
        assert response.status == 200
        // The analysis result as a set of Map objects.
        def actual = response.data.collect {
            [pathway: it.pathway, regulationType: it.regulationType]
        } as Set
        // Compare the matching set members. There might be new pathways
        // added over time and the scores might change, but the
        // regulation type should be the same. Therefore, validate only
        // that the expected pathways are in the actual result and
        // conserve the regulation type.
        assert actual.containsAll(expected)
        // Assertion passed.
        return true
    }
}
