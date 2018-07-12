package org.reactome.gsea.service;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URISyntaxException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.commons.collections4.IteratorUtils;
import org.apache.commons.collections4.Transformer;
import org.apache.commons.collections4.iterators.TransformIterator;
import org.reactome.gsea.config.PreRanked;
import org.reactome.gsea.model.GseaAnalysisResult;
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
import edu.mit.broad.genome.objects.GeneSet;
import edu.mit.broad.genome.objects.GeneSetMatrix;
import edu.mit.broad.genome.objects.RankedList;
import edu.mit.broad.genome.objects.esmatrix.db.EnrichmentDb;
import edu.mit.broad.genome.objects.esmatrix.db.EnrichmentResult;
import edu.mit.broad.genome.objects.esmatrix.db.EnrichmentScore;
import edu.mit.broad.genome.parsers.GmtParser;
import xtools.api.AbstractTool.Helper;
import xtools.api.param.GeneSetScoringTableReqdParam;
import xtools.api.param.IntegerParam;
import xtools.api.param.NormModeReqdParam;
import xtools.api.param.ParamFactory;

/**
 * @author Fred Loney <loneyf@ohsu.edu>
 */
// GSEA does not parameterize types, so we compensate with unchecked casts.
@SuppressWarnings("unchecked")
@RestController
public class GseaController {

    private static final String GMT_RESOURCE = "ReactomePathways.gmt";
    private static final String REACTOME_GMT = "ReactomePathways";

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
            @RequestParam(value="nperms", required=false) Integer nperms,
            @RequestParam(value="dataSetSizeMin", required=false) Integer dataSetSizeMin,
            @RequestParam(value="dataSetSizeMax", required=false) Integer dataSetSizeMax,
            @RequestBody List<List<String>> payload
    )
            throws Exception {
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
        return analyseRanked(nperms, dataSetSizeMin, dataSetSizeMax, rankedList);
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
    public @ResponseBody List<GseaAnalysisResult> analyseText(
            @RequestParam(value="nperms", required=false) Integer nperms,
            @RequestParam(value="dataSetSizeMin", required=false) Integer dataSetSizeMin,
            @RequestParam(value="dataSetSizeMax", required=false) Integer dataSetSizeMax,
            @RequestBody String payload
    )
            throws Exception {
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
 
        return analyseRanked(nperms, dataSetSizeMin, dataSetSizeMax, rankedList);
    }

    private List<GseaAnalysisResult> analyseRanked(Integer nperms,
            Integer dataSetSizeMin, 
            Integer dataSetSizeMax, 
            RankedList rankedList)
            throws Exception {
        GeneSet[] geneSets = getGeneSets(rankedList, dataSetSizeMin, dataSetSizeMax);
        EnrichmentResult[] erArray = analyse(rankedList, geneSets, nperms);
        List<EnrichmentResult> ers = Arrays.asList(erArray);
        // Work around the following bug:
        // * GSEA geneset parser converts the pathway names to upper case.
        // The work-around is to reread the gmt file into an {upper: lower}
        // map and replace the parsed GSEA gene set names.
        Set<String> erNames = new HashSet<String>();
        ers.stream().forEach(er -> erNames.add(er.getGeneSetName()));
        // The upper: lower map.
        final Map<String, String> utol = new HashMap<String, String>(erArray.length);
        // The lower: stable id map.
        final Map<String, String> stableIdMap = new HashMap<String, String>(erArray.length);
        ers.stream().forEach(er -> erNames.add(er.getGeneSetName()));
        
        initPathwayNameMaps(erNames, utol, stableIdMap);
 
        Transformer<EnrichmentResult, GseaAnalysisResult> erXfm =
                new Transformer<EnrichmentResult, GseaAnalysisResult>() {
            @Override
            public GseaAnalysisResult transform(EnrichmentResult er) {
                GseaAnalysisResult result = new GseaAnalysisResult();
                String ucName = er.getGeneSetName();
                GseaAnalysisResult.Pathway pathway = new GseaAnalysisResult.Pathway();
                String lcName = utol.get(ucName);
                pathway.setName(lcName);
                pathway.setStId(stableIdMap.get(lcName));
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
    
    private void initPathwayNameMaps(Set<String> erNames,
            Map<String, String> upCaseToLowCase,
            Map<String, String> nameToStableId) throws IOException {
        InputStream resource = getClass().getClassLoader().getResourceAsStream(GMT_RESOURCE);
        InputStreamReader reader = new InputStreamReader(resource);
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line = br.readLine()) != null) {
            String[] fields = line.split("\t");
            String lc = fields[0];
            String uc = lc.toUpperCase();
            if (erNames.contains(uc)) {
                upCaseToLowCase.put(uc, lc);
                String stId = fields[1];
                nameToStableId.put(lc, stId);
            }
        }
        br.close();
        reader.close();
        resource.close();
    }

    private RankedList getRankedList(String payload) {
        // Parse the input.
        String[] lines = payload.split(System.getProperty("line.separator"));
        List<List<String>> parsed = Stream.of(lines)
                                           .map(line -> Arrays.asList(line.split("\t")))
                                           .collect(Collectors.toList());
        return getRankedList(parsed);
    }

    private EnrichmentResult[] analyse(
            RankedList rankedList, GeneSet[] geneSets, Integer npermsOpt) throws Exception {
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

    private GeneSet[] getGeneSets(RankedList rankedList,
                                  Integer dataSetSizeMinOpt,
                                  Integer dataSetSizeMaxOpt) throws Exception {
        GmtParser gmtParser = new GmtParser();
//        InputStream gmtis = new FileInputStream(GMT_PATH);
        InputStream gmtis = getClass().getClassLoader().getResourceAsStream(GMT_RESOURCE);
        // The GSEA .gmt file parser returns a singleton array
        // consisting of the geneset matrix.
        List<GeneSetMatrix> gsms = (List<GeneSetMatrix>)gmtParser.parse(REACTOME_GMT, gmtis);
        GeneSetMatrix gsm = gsms.get(0);
        // The parsed gene sets.
        GeneSet[] gsmGeneSets = gsm.getGeneSets();
        int dataSetSizeMin = dataSetSizeMinOpt == null ? PreRanked.DEF_MIN_DATASET_SIZE : dataSetSizeMinOpt.intValue();
        int dataSetSizeMax = dataSetSizeMaxOpt == null ? PreRanked.DEF_MAX_DATASET_SIZE : dataSetSizeMaxOpt.intValue();
        IntegerParam fGeneSetMinSizeParam = ParamFactory.createGeneSetMinSizeParam(dataSetSizeMin, false);
        IntegerParam fGeneSetMaxSizeParam = ParamFactory.createGeneSetMaxSizeParam(dataSetSizeMax, false);
        // GSEA routine to filter the gene sets based on the min/max parameters.
        GeneSet[] geneSets = Helper.getGeneSets(rankedList, gsmGeneSets, fGeneSetMinSizeParam , fGeneSetMaxSizeParam);

        return geneSets;
    }
}
