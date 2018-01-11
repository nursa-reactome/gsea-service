package org.reactome.gsea.service;


import java.util.Arrays;
import java.util.Map.Entry;
import java.util.stream.Collectors;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
 
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Component;
import org.springframework.web.servlet.HandlerInterceptor;
import org.springframework.web.servlet.ModelAndView;
 
@Component
public class LogInterceptor implements HandlerInterceptor {
 
    Logger log = LoggerFactory.getLogger(this.getClass());
 
    @Override
    public void postHandle(HttpServletRequest request, HttpServletResponse response, Object object, ModelAndView model)
            throws Exception {
        String path = request.getServletPath().substring(1);
        log.info("Request \"" + path + "\" completed with status " + response.getStatus() + ".");
    }
 
    @Override
    public boolean preHandle(HttpServletRequest request, HttpServletResponse response, Object object) throws Exception {
        String params = request.getParameterMap().entrySet().stream()
                                 .map(LogInterceptor::formatParameter)
                                 .collect(Collectors.joining(", "));
        String path = request.getServletPath().substring(1);
        log.info("Issuing request \"" + path + "\" with parameters: {" + params + "}...");
        return true;
    }
    
    private static String formatParameter(Entry<String, String[]> entry) {
        String[] value = entry.getValue();
        String valsStr = Arrays.asList(value).stream()
                            .map(s -> "\"" + s + "\"")
                            .collect(Collectors.joining(","));
        String valueStr = value.length == 1 ? valsStr : "[" + valsStr + "]";
        
        return entry.getKey() + ": " + valueStr;
    }
 
}
