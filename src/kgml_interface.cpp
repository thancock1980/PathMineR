#include <libxml/xmlreader.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>

#include <string.h>
#include <iostream>
using namespace std;


#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>



extern "C" SEXP readkgmlfile(SEXP FILENAME) {
	const char *filename = CHAR(STRING_ELT(FILENAME,0));

	xmlDocPtr doc;
	xmlXPathContextPtr xpathCtx;
	xmlXPathObjectPtr nodes;

	//Check if the xml file has a KEGG DTD System.

	/* Load XML document */
	doc = xmlParseFile(filename);
	if (doc == NULL) {
		fprintf(stderr, "Error: unable to parse file \"%s\"\n", filename);
		return(R_NilValue);
	}

	/* Create xpath evaluation context */
	xpathCtx = xmlXPathNewContext(doc);
	if(xpathCtx == NULL) {
		fprintf(stderr,"Error: unable to create new XPath context\n");
		xmlFreeDoc(doc);
		return(R_NilValue);
	}

	/* Evaluate xpath expression */
	nodes = xmlXPathEvalExpression((xmlChar *) "//reaction", xpathCtx);
	if(nodes == NULL) {
		fprintf(stderr,"Error: unable to evaluate xpath expression \"%s\"\n", "//reaction");
		xmlXPathFreeContext(xpathCtx);
		xmlFreeDoc(doc);
		return(R_NilValue);
	}


	/* Parse XML Reactions*/
	xmlNodePtr curReaction;
    int size;
    int i;
    size = (nodes->nodesetval) ? nodes->nodesetval->nodeNr : 0;

    SEXP REACTIONLIST,ID;
    PROTECT(REACTIONLIST = allocVector(VECSXP,size));
    PROTECT(ID = allocVector(STRSXP,size));


    for(i = 0; i < size; ++i) {
    	curReaction = nodes->nodesetval->nodeTab[i];
    	xpathCtx->node = curReaction;

    	const char *name = (char *) xmlGetProp(curReaction,(const xmlChar *)"name");
    	SET_STRING_ELT(ID,i,mkChar(name));

    	SEXP REACTION,REACTIONNAMES;
		SEXP NAME, REVERSIBLE, REACTANTS, RSTOIC, PRODUCTS, PSTOIC, GENES, KINETICS,KNAMES;

		PROTECT(REACTION = allocVector(VECSXP,8));
		PROTECT(REACTIONNAMES = allocVector(STRSXP,8));


		PROTECT(NAME = allocVector(STRSXP,1));
		SET_STRING_ELT(NAME,0, mkChar(name));
		SET_VECTOR_ELT(REACTION,0,NAME); SET_STRING_ELT(REACTIONNAMES,0,mkChar("name"));

		PROTECT(REVERSIBLE = allocVector(LGLSXP,1));
		cout << (char *) xmlGetProp(curReaction,(const xmlChar *)"type") <<endl;
		LOGICAL(REVERSIBLE)[0] = strcmp((char *) xmlGetProp(curReaction,(const xmlChar *)"type"), "irreversible");
		SET_VECTOR_ELT(REACTION,1,REVERSIBLE); SET_STRING_ELT(REACTIONNAMES,1,mkChar("reversible"));


		xmlNodeSetPtr reactantNodes = xmlXPathEvalExpression( (const xmlChar *) "./substrate", xpathCtx)->nodesetval;
		int numOfReactants = (reactantNodes) ? reactantNodes->nodeNr : 0;
		PROTECT(REACTANTS = allocVector(STRSXP,numOfReactants));
		PROTECT(RSTOIC = allocVector(REALSXP,numOfReactants));

		for (int r = 0;r < numOfReactants;r++) {
			SET_STRING_ELT(REACTANTS,r,mkChar((char*) xmlGetProp(reactantNodes->nodeTab[r],(const xmlChar *)"name")));
			REAL(RSTOIC)[r] = NA_REAL;
		}
		SET_VECTOR_ELT(REACTION,2,REACTANTS); SET_STRING_ELT(REACTIONNAMES,2,mkChar("reactants"));
		SET_VECTOR_ELT(REACTION,3,RSTOIC);SET_STRING_ELT(REACTIONNAMES,3,mkChar("reactant.stoichiometry"));
		cout << "Reactant level :|" << endl;

		xmlNodeSetPtr productNodes = xmlXPathEvalExpression( (const xmlChar *) "./product", xpathCtx)->nodesetval;
		int numOfProducts = (productNodes) ? productNodes->nodeNr : 0;
		PROTECT(PRODUCTS = allocVector(STRSXP,numOfProducts));
		PROTECT(PSTOIC = allocVector(REALSXP,numOfProducts));

		for (int p = 0;p < numOfProducts;p++) {
			SET_STRING_ELT(PRODUCTS,p,mkChar((char*) xmlGetProp(productNodes->nodeTab[p],(const xmlChar *)"name")));
			REAL(PSTOIC)[p] = NA_REAL;
		}
		SET_VECTOR_ELT(REACTION,4,PRODUCTS); SET_STRING_ELT(REACTIONNAMES,4,mkChar("products"));
		SET_VECTOR_ELT(REACTION,5,PSTOIC);SET_STRING_ELT(REACTIONNAMES,5,mkChar("product.stoichiometry"));
		cout << "Product level :|" << endl;

        SET_VECTOR_ELT(REACTION,6,R_NilValue);SET_STRING_ELT(REACTIONNAMES,6,mkChar("kinetics"));


        string geneXPath = ((string)"//entry[@type='gene' and @reaction='")+((string)name)+((string)"']");
        xmlNodeSetPtr geneNodes = xmlXPathEvalExpression(
        		(const xmlChar *) geneXPath.c_str(), xpathCtx)->nodesetval;
        int numOfModifiers = (geneNodes) ? geneNodes->nodeNr : 0;
		PROTECT(GENES = allocVector(STRSXP,numOfModifiers));
		for (int m = 0;m < numOfModifiers;m++) {
			SET_STRING_ELT(GENES,m,mkChar((char*) xmlGetProp(geneNodes->nodeTab[m],(const xmlChar *)"name")));
		}
		SET_VECTOR_ELT(REACTION,7,GENES); SET_STRING_ELT(REACTIONNAMES,7,mkChar("genes"));

		setAttrib(REACTION,R_NamesSymbol,REACTIONNAMES);
        SET_VECTOR_ELT(REACTIONLIST,i,REACTION);

        UNPROTECT(9);

    }

    setAttrib(REACTIONLIST,R_NamesSymbol,ID);
    UNPROTECT(2);

	/* Cleanup */
	xmlXPathFreeObject(nodes);
	xmlXPathFreeContext(xpathCtx);
	xmlFreeDoc(doc);

    //cout << "RacList returns :|" << endl;
    return(REACTIONLIST);
};

/*
extern "C" SEXP readkgmlfile(SEXP FILENAME, SEXP XPATH) {
	const char *filename = CHAR(STRING_ELT(FILENAME,0));
	const char *x = CHAR(STRING_ELT(XPATH,0));
	execute_xpath_expression(filename, BAD_CAST x);

	return(R_NilValue);
}


void parseTree(const char *filename){
	xmlDoc *doc = NULL;
	xmlNode *root_element = NULL;
	xmlNode *cur_node = NULL;

	doc = xmlReadFile(filename, NULL, 0);
	root_element = xmlDocGetRootElement(doc);

	for (cur_node = root_element->children; cur_node; cur_node = cur_node->next) {
			printf("node type: Element, name: %s\n", cur_node->name, cur_node->type);
	}

}

void streamFile(const char *filename) {
    xmlTextReaderPtr reader;
    int ret;

    reader = xmlReaderForFile(filename, NULL, 0);
    if (reader != NULL) {
        ret = xmlTextReaderRead(reader);
        while (ret == 1) {
            processNode(reader);
            ret = xmlTextReaderRead(reader);
        }
        xmlFreeTextReader(reader);
        if (ret != 0) {
            fprintf(stderr, "%s : failed to parse\n", filename);
        }
    } else {
        fprintf(stderr, "Unable to open %s\n", filename);
    }
}

void processNode(xmlTextReaderPtr reader) {
    const xmlChar *name, *value;

    name = xmlTextReaderConstName(reader);
    if (name == NULL)
	name = BAD_CAST "--";

    value = xmlTextReaderConstValue(reader);

    printf("%d %d %s %d %d",
	    xmlTextReaderDepth(reader),
	    xmlTextReaderNodeType(reader),
	    name,
	    xmlTextReaderIsEmptyElement(reader),
	    xmlTextReaderHasValue(reader));
    if (value == NULL)
	printf("\n");
    else {
        if (xmlStrlen(value) > 40)
            printf(" %.40s...\n", value);
        else
	    printf(" %s\n", value);
    }
}




/*
void
print_xpath_nodes(xmlNodeSetPtr nodes, FILE* output) {
    xmlNodePtr cur;
    int size;
    int i;


    size = (nodes) ? nodes->nodeNr : 0;

    fprintf(output, "Result (%d nodes):\n", size);
    for(i = 0; i < size; ++i) {


	if(nodes->nodeTab[i]->type == XML_NAMESPACE_DECL) {
	    xmlNsPtr ns;

	    ns = (xmlNsPtr)nodes->nodeTab[i];
	    cur = (xmlNodePtr)ns->next;
	    if(cur->ns) {
	        fprintf(output, "= namespace \"%s\"=\"%s\" for node %s:%s\n",
		    ns->prefix, ns->href, cur->ns->href, cur->name);
	    } else {
	        fprintf(output, "= namespace \"%s\"=\"%s\" for node %s\n",
		    ns->prefix, ns->href, cur->name);
	    }
	} else if(nodes->nodeTab[i]->type == XML_ELEMENT_NODE) {
	    cur = nodes->nodeTab[i];
	    if(cur->ns) {
    	        fprintf(output, "= element node \"%s:%s\"\n",
		    cur->ns->href, cur->name);
	    } else {
    	        fprintf(output, "= element node \"%s\"\n",
		    cur->name);
	    }
	} else {
	    cur = nodes->nodeTab[i];
	    fprintf(output, "= node \"%s\": type %d\n", cur->name, cur->type);
	}
    }
}
*/
