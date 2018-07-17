package org.reactome.gsea.service;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.httpclient.HttpClient;
import org.apache.commons.httpclient.HttpMethod;
import org.apache.commons.httpclient.HttpStatus;
import org.apache.commons.httpclient.methods.GetMethod;
import org.apache.commons.httpclient.methods.PostMethod;
import org.apache.commons.httpclient.methods.RequestEntity;
import org.apache.commons.httpclient.methods.StringRequestEntity;
import org.junit.Test;

import com.fasterxml.jackson.databind.ObjectMapper;

public class WSTest {
    protected final String HOST_URL = "http://localhost:8070/";
//    private final String HOST_URL = "http://reactomews.oicr.on.ca:8080/gsea-service/";
    private final String HTTP_POST = "Post";
    private final String HTTP_GET = "Get";
    
    public WSTest() {
    }
    
    @Test
    public void testAnalysis() throws Exception {
        String url = HOST_URL + "analyse";
        InputStream resource = this.getClass().getResourceAsStream("/fixtures/gTqItVnDEP.rnk");
        InputStreamReader reader = new InputStreamReader(resource);
        BufferedReader br = new BufferedReader(reader);
        List<List<String>> lines = new ArrayList<>();
        String line = null;
        while ((line = br.readLine()) != null) {
            String[] tokens = line.split("(\t| )+");
            lines.add(Arrays.asList(tokens));
        }
        br.close();
        reader.close();
        resource.close();
        System.out.println();
        
        String query = convertListToString(lines);
        
        String results = callHttp(url, HTTP_POST, "text/plain", query);
        System.out.println("Analysis results via query in text/plain:\n" + results);
        
        ObjectMapper mapper = new ObjectMapper();
        query = mapper.writeValueAsString(lines);
        results = callHttp(url, HTTP_POST, "application/json", query);
        System.out.println("\nAnalysis results via query in application/json:\n" + results);

    }
    
    private String convertListToString(List<List<String>> lines) {
        StringBuilder builder = new StringBuilder();
        lines.forEach(line -> {
            line.forEach(token -> builder.append(token).append("\t"));
            builder.append("\n");
        });
        return builder.toString();
    }
    
    protected String callHttp(String url,
            String type,
            String requestType,
            String query) throws IOException {
        HttpMethod method = null;
        HttpClient client = null;
        if (type.equals(HTTP_POST)) {
            method = new PostMethod(url);
            client = initializeHTTPClient((PostMethod) method, requestType, query);
        } else {
            method = new GetMethod(url); // Default
            client = new HttpClient();
        }
        method.setRequestHeader("Accept", "application/json");
        int responseCode = client.executeMethod(method);
        if (responseCode == HttpStatus.SC_OK) {
            InputStream is = method.getResponseBodyAsStream();
            return readMethodReturn(is);
        } else {
            System.err.println("Error from server: " + method.getResponseBodyAsString());
            System.out.println("Response code: " + responseCode);
            throw new IllegalStateException(method.getResponseBodyAsString());
        }
    }
    
    protected String readMethodReturn(InputStream is) throws IOException {
        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader reader = new BufferedReader(isr);
        StringBuilder builder = new StringBuilder();
        String line = null;
        while ((line = reader.readLine()) != null)
            builder.append(line).append("\n");
        reader.close();
        isr.close();
        is.close();
        // Remove the last new line
        String rtn = builder.toString();
        // Just in case an empty string is returned
        if (rtn.length() == 0)
            return rtn;
        return rtn.substring(0, rtn.length() - 1);
    }
    
    private HttpClient initializeHTTPClient(PostMethod post, 
                                            String requestType,
                                            String query) throws UnsupportedEncodingException {
        RequestEntity entity = new StringRequestEntity(query, requestType, "UTF-8");
//        RequestEntity entity = new StringRequestEntity(query, "application/json", "UTF-8");
        post.setRequestEntity(entity);
//        post.setRequestHeader("Accept", "application/JSON, application/XML, text/plain");
              post.setRequestHeader("Accept", "application/json");
        HttpClient client = new HttpClient();
        return client;
    }

}
