package org.reactome.gsea.controller;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.commons.collections4.IteratorUtils;
import org.apache.commons.collections4.Transformer;
import org.apache.commons.collections4.iterators.TransformIterator;
import org.apache.log4j.Logger;
import org.reactome.gsea.config.PreRanked;
import org.reactome.gsea.model.GseaAnalysisResult;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.context.annotation.PropertySource;
import org.springframework.core.env.Environment;
import org.springframework.web.bind.annotation.RequestBody;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestMethod;
import org.springframework.web.bind.annotation.RequestParam;
import org.springframework.web.bind.annotation.ResponseBody;
import org.springframework.web.bind.annotation.RestController;

import edu.mit.broad.genome.alg.RankedListGenerators;
import edu.mit.broad.genome.alg.gsea.GeneSetCohortGenerator;
import edu.mit.broad.genome.alg.gsea.KSTests;
import edu.mit.broad.genome.alg.gsea.PValueCalculator;
import edu.mit.broad.genome.alg.gsea.PValueCalculatorImpls;
import edu.mit.broad.genome.math.Order;
import edu.mit.broad.genome.math.RandomSeedGenerator;
import edu.mit.broad.genome.math.RandomSeedGenerators;
import edu.mit.broad.genome.math.SortMode;
import edu.mit.broad.genome.objects.FSet;
import edu.mit.broad.genome.objects.GeneSet;
import edu.mit.broad.genome.objects.RankedList;
import edu.mit.broad.genome.objects.esmatrix.db.EnrichmentDb;
import edu.mit.broad.genome.objects.esmatrix.db.EnrichmentResult;
import edu.mit.broad.genome.objects.esmatrix.db.EnrichmentScore;
import xtools.api.AbstractTool.Helper;
import xtools.api.param.GeneSetScoringTableReqdParam;
import xtools.api.param.IntegerParam;
import xtools.api.param.NormModeReqdParam;
import xtools.api.param.ParamFactory;

/**
 * @author Fred Loney <loneyf@ohsu.edu>
 */
// GSEA does not parameterize types, so we compensate with unchecked casts.
@PropertySource("classpath:application.properties")
@RestController
public class GseaController {
    private static final Logger logger = Logger.getLogger(GseaController.class);
    @Autowired
    private Environment env;
    // Cached loaded GeneSet[] to save some time 
    private GeneSet[] humanGeneSets;
    private GeneSet[] mouseGeneSets;
    
    private String getGmtResource(String resourceName) {
        String gmtResource = env.getProperty(resourceName);
        if (gmtResource == null)
            throw new IllegalStateException("gmtResource has not been set!");
        logger.debug("Set gmtResource as: " + gmtResource);
        return gmtResource;
    }
    
    /**
     * Runs a GSEA analysis on a preranked list.
     * 
     * @param payload the {gene symbol, rank value} JSON records
     *  request body
     * @return the GSEA execution result response body
     * @throws URISyntaxException
     * @throws IOException
     */
    @RequestMapping(value="/analyse", 
                    method=RequestMethod.POST,
                    consumes = "application/json")
    public @ResponseBody List<GseaAnalysisResult> analyse(
            @RequestParam(value="species", required=false, defaultValue = "human") String species,
            @RequestParam(value="nperms", required=false, defaultValue = "1000") Integer nperms,
            @RequestParam(value="dataSetSizeMin", required=false, defaultValue = "3") Integer dataSetSizeMin,
            @RequestParam(value="dataSetSizeMax", required=false, defaultValue = "500") Integer dataSetSizeMax,
            @RequestBody List<List<String>> payload) throws Exception {
        // This method and the methods it calls are adapted from the GSEA
        // Preranked internal implementation.
        // There is no clean and simple GSEA public API for performing the
        // enrichment on this method's parameters. Rather, there is an implicit
        // GSEA assumption that Preranked is called with input file name parameters
        // and produces an output report. GSEA is deeply entwined with these
        // assumptions. Bits and pieces in the methods below are borrowed from
        // different places in the GSEA code to replicate the GSEA file parsing
        // and analysis execution without file parsing and report production.
        RankedList rankedList = getRankedList(payload);
        return analyseRanked(species, nperms, dataSetSizeMin, dataSetSizeMax, rankedList);
    }

    /**
     * Runs a GSEA analysis on a preranked list.
     * 
     * @param payload the {gene symbol, rank value} JSON records
     *  request body
     * @return the GSEA execution result response body
     * @throws URISyntaxException
     * @throws IOException
     */
    @RequestMapping(value="/analyse",
                    method=RequestMethod.POST,
                    consumes = "text/plain")
    public @ResponseBody List<GseaAnalysisResult> analyse(
            @RequestParam(value="species", required=false, defaultValue = "human") String species,
            @RequestParam(value="nperms", required=false, defaultValue = "1000") Integer nperms,
            @RequestParam(value="dataSetSizeMin", required=false, defaultValue = "3") Integer dataSetSizeMin,
            @RequestParam(value="dataSetSizeMax", required=false, defaultValue = "500") Integer dataSetSizeMax,
            @RequestBody String payload) throws Exception {
        RankedList rankedList = getRankedList(payload);
 
        return analyseRanked(species, nperms, dataSetSizeMin, dataSetSizeMax, rankedList);
    }

    private List<GseaAnalysisResult> analyseRanked(String species,
                                                   Integer nperms,
                                                   Integer dataSetSizeMin, 
                                                   Integer dataSetSizeMax, 
                                                   RankedList rankedList) throws Exception {
        GeneSet[] geneSets = getGeneSets(species, rankedList, dataSetSizeMin, dataSetSizeMax);
        EnrichmentResult[] erArray = analyse(rankedList, geneSets, nperms);
        List<EnrichmentResult> ers = Arrays.asList(erArray);
        // Convert into a simple result for json
        Transformer<EnrichmentResult, GseaAnalysisResult> erXfm =
                new Transformer<EnrichmentResult, GseaAnalysisResult>() {
            @Override
            public GseaAnalysisResult transform(EnrichmentResult er) {
                GseaAnalysisResult result = new GseaAnalysisResult();
                GeneSet geneset = er.getGeneSet();
                GseaAnalysisResult.Pathway pathway = new GseaAnalysisResult.Pathway();
                pathway.setName(geneset.getName());
                pathway.setStId(geneset.getNameEnglish()); // This is stable id (the second column. See parseGeneSet).
                result.setPathway(pathway);
                EnrichmentScore es = er.getScore();
                result.setHitCount(es.getNumHits());
                result.setScore(es.getES());
                result.setNormalizedScore(es.getNES());
                result.setPvalue(es.getNP());
                result.setFdr(es.getFDR());
                return result;
            }
        };
        Iterator<GseaAnalysisResult> erIter =
                new TransformIterator<EnrichmentResult, GseaAnalysisResult>(ers.iterator(), erXfm);
        List<GseaAnalysisResult> result = IteratorUtils.toList(erIter);
 
        return result;
    }

    private RankedList getRankedList(String payload) {
        // Parse the input.
        String[] lines = payload.split(System.getProperty("line.separator"));
        List<List<String>> parsed = Stream.of(lines)
                                           .map(line -> Arrays.asList(line.split("\t")))
                                           .collect(Collectors.toList());
        return getRankedList(parsed);
    }

    private EnrichmentResult[] analyse(RankedList rankedList, 
                                       GeneSet[] geneSets,
                                       Integer npermsOpt) throws Exception {
        RandomSeedGenerator rst = RandomSeedGenerators.lookup("timestamp");
        // The cohort generator is filled in and used internally by GSEA.
        GeneSetScoringTableReqdParam fGcohGenReqdParam = new GeneSetScoringTableReqdParam();
        GeneSetCohortGenerator gcohgen = fGcohGenReqdParam.createGeneSetCohortGenerator(true);
        // The default number of permutations is 1000.
        int nperms = npermsOpt == null ? 1000 : npermsOpt.intValue();
        // The GSEA internal KS test implementation.
        KSTests tests = new KSTests();
        // There are two enrichment phases. The first phase calculates the base score.
        // The second phase calculates the normalized score, FDR and p-values.
        EnrichmentDb edbPhaseOne = tests.executeGsea(
                    rankedList,
                    geneSets,
                    nperms, // permutations
                    rst,
                    null, // unused chip argument
                    gcohgen
            );
        
        // The GSEA normalization mode parameter.
        NormModeReqdParam modeParam = new NormModeReqdParam();
        String modeDef = (String) modeParam.getDefault();
        // Phase two.
        PValueCalculator pvc = new PValueCalculatorImpls.GseaImpl(modeDef);
        EnrichmentResult[] results = pvc.calcNPValuesAndFDR(edbPhaseOne.getResults());

        return results;
    }

    private RankedList getRankedList(List<List<String>> payload) {
        // Use a map to avoid any duplication in the payload
        Map<String, Float> nameToValue = new HashMap<>();
        payload.stream().forEach(pair -> {
            nameToValue.compute(pair.get(0), (key, value) -> {
                return new Float(pair.get(1));
            });
        });
        String[] names =
                nameToValue.keySet().toArray(new String[nameToValue.size()]);
        float[] values = new float[nameToValue.size()];
        for (int i=0; i < names.length; i++) {
            values[i] = nameToValue.get(names[i]);
        }
        
        return RankedListGenerators.createBySorting(
                "REST", names, values, SortMode.REAL, Order.DESCENDING);
    }

    private GeneSet[] getGeneSets(String species,
                                  RankedList rankedList,
                                  Integer dataSetSizeMinOpt,
                                  Integer dataSetSizeMaxOpt) throws Exception {
        GeneSet[] gsmGeneSets = getGeneSets(species);
        int dataSetSizeMin = dataSetSizeMinOpt == null ? PreRanked.DEF_MIN_DATASET_SIZE : dataSetSizeMinOpt.intValue();
        int dataSetSizeMax = dataSetSizeMaxOpt == null ? PreRanked.DEF_MAX_DATASET_SIZE : dataSetSizeMaxOpt.intValue();
        IntegerParam fGeneSetMinSizeParam = ParamFactory.createGeneSetMinSizeParam(dataSetSizeMin, false);
        IntegerParam fGeneSetMaxSizeParam = ParamFactory.createGeneSetMaxSizeParam(dataSetSizeMax, false);
        // GSEA routine to filter the gene sets based on the min/max parameters.
        GeneSet[] geneSets = Helper.getGeneSets(rankedList, gsmGeneSets, fGeneSetMinSizeParam , fGeneSetMaxSizeParam);

        return geneSets;
    }

    private GeneSet[] getGeneSets(String species) throws Exception {
        if (species.equals("mouse"))
            return getMouseGeneSets();
        return getHumanGeneSets(); // As the default
    }
    
    private GeneSet[] getMouseGeneSets() throws Exception {
        if (mouseGeneSets != null)
            return mouseGeneSets;
        String resource = getGmtResource("mouseGmtResource");
        mouseGeneSets = parseGmtFile(resource);
        return mouseGeneSets;
    }
    
    private GeneSet[] getHumanGeneSets() throws Exception {
        if (humanGeneSets != null)
            return humanGeneSets;
        String resource = getGmtResource("gmtResource");
        humanGeneSets = parseGmtFile(resource);
        return humanGeneSets;
    }
    
    /**
     * Modified from GmtParser to have a lean parser without upper case pathway names.
     * @param resource
     * @return
     * @throws IOException
     */
    private GeneSet[] parseGmtFile(String resource) throws IOException {
        logger.debug("Loading gene sets for " + resource);
        InputStream gmtis = getClass().getClassLoader().getResourceAsStream(resource);
        BufferedReader br = new BufferedReader(new InputStreamReader(gmtis));
        String line = null;
        List<GeneSet> genesets = new ArrayList<>();
        while ((line = br.readLine()) != null) {
            String[] tokens = line.split("\t");
            List<String> geneNames = Arrays.asList(tokens)
                                           .subList(2, tokens.length);
            GeneSet gset = new FSet(tokens[0],
                                    tokens[1],
                                    geneNames,
                                    true);
            genesets.add(gset);
        }
        br.close();
        logger.debug("Total gene sets: " + genesets.size());
        return genesets.toArray(new GeneSet[] {});
    }
    
}
