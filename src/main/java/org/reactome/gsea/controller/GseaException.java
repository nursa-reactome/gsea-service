package org.reactome.gsea.controller;

public class GseaException extends RuntimeException {

    private static final long serialVersionUID = -5346430602050386666L;

    public GseaException(String string) {
        super(string);
    }

    public GseaException(String string, Exception e) {
        super(string, e);
    }

}
