package org.reactome.gsea.service;

public final class GseaException extends RuntimeException {

    private static final long serialVersionUID = 5785726530598739074L;

    public GseaException(String message, Exception cause) {
        super(message, cause);
    }
    
}
