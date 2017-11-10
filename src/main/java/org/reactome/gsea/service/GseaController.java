package org.reactome.gsea.service;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.net.URISyntaxException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.collections4.IteratorUtils;
import org.apache.commons.collections4.Transformer;
import org.apache.commons.collections4.iterators.TransformIterator;
import org.reactome.gsea.model.AnalysisResult;
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

    private static final String GMT_PATH = "/usr/local/reactomes/Reactome/production/GKB/scripts/release/WebELVTool/ReactomePathways.gmt";

    /**
     * Runs a GSEA analysis on a preranked list.
     * 
     * @param rankedList the tab-separated {gene symbol, rank value} records
     *  request body
     * @return the GSEA execution result response body
     * @throws URISyntaxException
     * @throws IOException
     */
    @RequestMapping(value = "/analyse", method = RequestMethod.POST)
    public @ResponseBody List<AnalysisResult> analyse(
           @RequestParam(value="nperms", required=false) Integer nperms,
            @RequestBody List<List<String>> payload
    )
            throws URISyntaxException, IOException {

        Transformer<EnrichmentResult, AnalysisResult> erXfm =
                new Transformer<EnrichmentResult, AnalysisResult>() {
            @Override
            public AnalysisResult transform(EnrichmentResult er) {
                AnalysisResult result = new AnalysisResult();
                result.setPathway(er.getGeneSetName());
                EnrichmentScore es = er.getScore();
                result.setScore(es.getES());
                result.setNormalizedScore(es.getNES());
                result.setPvalue(es.getNP());
                result.setFdr(es.getFDR());
                return result;
            }
        };

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
        GeneSet[] geneSets = getGeneSets(rankedList);
        EnrichmentResult[] erArray = analyse(rankedList, geneSets, nperms);
        List<EnrichmentResult> ers = Arrays.asList(erArray);
        Iterator<AnalysisResult> erIter =
                new TransformIterator<EnrichmentResult, AnalysisResult>(ers.iterator(), erXfm);
        List<AnalysisResult> result = IteratorUtils.toList(erIter);
        List<String> resn =
                result.stream()
                      .map(AnalysisResult::getPathway)
                      .collect(Collectors.toList());
        System.out.println(">> n " + resn);
        List<Float> resv =
                result.stream()
                      .map(AnalysisResult::getScore)
                      .map(Float::new)
                      .collect(Collectors.toList());
        System.out.println(">> v " + resv);

        return result;
    }

    private EnrichmentResult[] analyse(RankedList rankedList, GeneSet[] geneSets, Integer npermsOpt) {
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
        EnrichmentDb edbPhaseOne;
        try {
            edbPhaseOne = tests.executeGsea(
                    rankedList,
                    geneSets,
                    nperms, // permutations
                    rst,
                    null, // unused chip argument
                    gcohgen
            );
        } catch (Exception e) {
            throw new GseaException("GSEA enrichment error", e);
        }
        // The GSEA normalization mode parameter.
        NormModeReqdParam modeParam = new NormModeReqdParam();
        String modeDef = (String) modeParam.getDefault();
        // Phase two.
        final PValueCalculator pvc = new PValueCalculatorImpls.GseaImpl(modeDef);
        final EnrichmentResult[] results = pvc.calcNPValuesAndFDR(edbPhaseOne.getResults());

        return results;
    }

    private RankedList getRankedList(List<List<String>> payload)
            throws UnsupportedEncodingException {
        // The current names/values index.
        int current = 0;
        // Name set to check duplicates.
        Set<String> nameSet = new HashSet<String>(payload.size());
        String[] names = new String[payload.size()];
        float[] values = new float[payload.size()];
        for (List<String> entry: payload) {
            String name = entry.get(0);
            if (!nameSet.contains(name)) {
                names[current] = name;
                String valueStr = entry.get(1);
                values[current] = Float.parseFloat(valueStr);
                nameSet.add(name);
                current++;
            }
        }
        // Resize the arrays if necessary to account for duplicate removal.
        if (current < payload.size()) {
            String[] namesUnique = new String[current];
            System.arraycopy(names, 0, namesUnique, 0, current);
            names = namesUnique;
            float[] valuesUnique = new float[current];
            System.arraycopy(values, 0, valuesUnique, 0, current);
            values = valuesUnique;
        }
        
        return RankedListGenerators.createBySorting(
                "REST", names, values, SortMode.REAL, Order.DESCENDING);
    }

    private GeneSet[] getGeneSets(RankedList rankedList) throws FileNotFoundException {
        GmtParser gmtParser = new GmtParser();
        FileInputStream gmtis = new FileInputStream(GMT_PATH);
        // The GSEA .gmt file parser returns a singleton array
        // consisting of the geneset matrix.
        List<GeneSetMatrix> gsms;
        try {
            gsms = (List<GeneSetMatrix>)gmtParser.parse(GMT_PATH, gmtis);
        } catch (Exception e) {
            throw new GseaException("GSEA gene lists parsing error", e);
        }
        GeneSetMatrix gsm = gsms.get(0);
        // The parsed gene sets.
        GeneSet[] gsmGeneSets = gsm.getGeneSets();
        // Default parameters lifted from AbstractGseaTool. 
        IntegerParam fGeneSetMinSizeParam = ParamFactory.createGeneSetMinSizeParam(15, false);
        IntegerParam fGeneSetMaxSizeParam = ParamFactory.createGeneSetMaxSizeParam(500, false);
        // GSEA routine to filter the gene sets based on the min/max parameters.
        GeneSet[] geneSets;
        try {
            geneSets = Helper.getGeneSets(rankedList, gsmGeneSets, fGeneSetMinSizeParam , fGeneSetMaxSizeParam);
        } catch (Exception e) {
            throw new GseaException("GSEA enrichment error", e);
        }

        return geneSets;
    }
}
