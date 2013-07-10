#include "sbml_interface.h"


extern "C" SEXP readsbmlfile(SEXP FILENAME) {
  SEXP SPECIESFRAME, REACTIONLIST, NETWORK, OUT,NAMES;
  const char *filename = CHAR(STRING_ELT(FILENAME,0));
    
  SBMLDocument* document = readSBML(filename);
  unsigned int errors = document->getNumErrors();
  
  cout << endl;
  cout << "Processing SBML file: :|" << filename << endl;
  cout << "\terror(s): " << errors  << endl;
  cout << "\tSBML level: " << document->getLevel() << endl;
  cout << "\tSBML version: " << document->getLevel() << endl;  
  if (errors > 0) document->printErrors(cerr);
  cout << endl;

  Model *model = document->getModel();
  SPECIESFRAME = getSpeciesFrame(model);
  REACTIONLIST = getReactionList(model);

  cout << "C++ returns :|" << endl;
  PROTECT(NETWORK = allocVector(STRSXP,1));
  SET_STRING_ELT(NETWORK,0,mkChar(model->getName().c_str()));  
  PROTECT(OUT = allocVector(VECSXP,3));
  PROTECT(NAMES = allocVector(STRSXP,3));
  SET_VECTOR_ELT(OUT,0,SPECIESFRAME); SET_STRING_ELT(NAMES,0,mkChar("species.info"));
  SET_VECTOR_ELT(OUT,1,REACTIONLIST); SET_STRING_ELT(NAMES,1,mkChar("reaction.info"));
  SET_VECTOR_ELT(OUT,2,NETWORK); SET_STRING_ELT(NAMES,2,mkChar("network"));
  setAttrib(OUT,R_NamesSymbol,NAMES);

  UNPROTECT(3);

  return(OUT);
};

SEXP getReactionList(Model *model) {
    ListOfReactions *reactions = model->getListOfReactions();
    
    cout << "Processing All Reactions" << endl;
    cout << "  Number of reactions found: " << reactions->size() << endl;
    
    SEXP REACTIONLIST,ID;
    PROTECT(REACTIONLIST = allocVector(VECSXP,reactions->size()));
    PROTECT(ID = allocVector(STRSXP,reactions->size()));
    

    for (unsigned i = 0;i < reactions->size();i++) {

    	//Get reaction stable ids
    	List cvl = List();
		RDFAnnotationParser::parseRDFAnnotation(reactions-> get(i)-> getAnnotation(), (List *) &cvl);
		CVTerm *cvt = (CVTerm *) cvl.get(0);
		unsigned int n=0;
		string URI;
		SEXP ST_ID;
		while(n< cvt -> getNumResources()){
			URI= cvt ->getResourceURI(n);
			int pos = URI.find("REACT_");
			if(pos > 0){
				PROTECT(ST_ID = mkChar(URI.substr(pos).c_str()));
				SET_STRING_ELT(ID,i, ST_ID);
				UNPROTECT(1);
				break;
			}
			n++;
		}

    	//SET_STRING_ELT(ID,i,mkChar(reactions->get(i)->getId().c_str()));

        SEXP REACTION,REACTIONNAMES;
        SEXP NAME, REVERSIBLE, REACTANTS, RSTOIC, PRODUCTS, PSTOIC, GENES, KINETICS,KNAMES;
    
        PROTECT(REACTION = allocVector(VECSXP,8));
        PROTECT(REACTIONNAMES = allocVector(STRSXP,8));
        
        PROTECT(NAME = allocVector(STRSXP,1));   
        SET_STRING_ELT(NAME,0,mkChar(reactions->get(i)->getName().c_str()));
        SET_VECTOR_ELT(REACTION,0,NAME); SET_STRING_ELT(REACTIONNAMES,0,mkChar("name"));
        
        PROTECT(REVERSIBLE = allocVector(LGLSXP,1));
        LOGICAL(REVERSIBLE)[0] = reactions->get(i)->getReversible();
        SET_VECTOR_ELT(REACTION,1,REVERSIBLE); SET_STRING_ELT(REACTIONNAMES,1,mkChar("reversible"));

        int numOfReactants = reactions->get(i)->getNumReactants();
        PROTECT(REACTANTS = allocVector(STRSXP,numOfReactants));  
        PROTECT(RSTOIC = allocVector(REALSXP,numOfReactants));
        for (int r = 0;r < numOfReactants;r++) {
            SET_STRING_ELT(REACTANTS,r,mkChar(reactions->get(i)->getReactant(r)->getSpecies().c_str()));
            REAL(RSTOIC)[r] = reactions->get(i)->getReactant(r)->getStoichiometry();
        }   
        SET_VECTOR_ELT(REACTION,2,REACTANTS); SET_STRING_ELT(REACTIONNAMES,2,mkChar("reactants"));
        SET_VECTOR_ELT(REACTION,3,RSTOIC);SET_STRING_ELT(REACTIONNAMES,3,mkChar("reactant.stoichiometry"));
        //cout << "Reactant level :|" << endl;
        
        int numOfProducts = reactions->get(i)->getNumProducts();
        PROTECT(PRODUCTS = allocVector(STRSXP,numOfProducts));   
        PROTECT(PSTOIC = allocVector(REALSXP,numOfProducts));
        for (int p = 0;p < numOfProducts;p++) {
            SET_STRING_ELT(PRODUCTS,p,mkChar(reactions->get(i)->getProduct(p)->getSpecies().c_str()));
            REAL(PSTOIC)[p] = reactions->get(i)->getProduct(p)->getStoichiometry();
        } 
        SET_VECTOR_ELT(REACTION,4,PRODUCTS); SET_STRING_ELT(REACTIONNAMES,4,mkChar("products"));
        SET_VECTOR_ELT(REACTION,5,PSTOIC);SET_STRING_ELT(REACTIONNAMES,5,mkChar("product.stoichiometry"));
        //cout << "Product level :|" << endl;
        //Kinetic law
        KineticLaw *kinetics = reactions->get(i)->getKineticLaw();
        int knum;
        if(kinetics!=NULL){knum = kinetics->getNumParameters();}
        else{knum = 0;}
		//cout << "Kinetic for level :|" << knum << endl;
		PROTECT(KINETICS = allocVector(VECSXP,knum));
        PROTECT(KNAMES = allocVector(STRSXP,knum));

        for (int k = 0;k < knum;k++) {
           SEXP value;
           PROTECT(value = allocVector(REALSXP,1));
           REAL(value)[0] = kinetics->getParameter(k)->getValue();
           SET_STRING_ELT(KNAMES,k,mkChar(kinetics->getParameter(k)->getId().c_str()));
           SET_VECTOR_ELT(KINETICS,k,value);
           UNPROTECT(1);
        }
        setAttrib(KINETICS,R_NamesSymbol,KNAMES);
        SET_VECTOR_ELT(REACTION,6,KINETICS);SET_STRING_ELT(REACTIONNAMES,6,mkChar("kinetics"));
        //cout << "kinetic level :|" << endl;

        // read genes from notes section of sbml
        /*string notes = reactions->get(i)->getNotesString();
        int g1 = notes.find("LOCUS:");
        int g2 = 0;
        vector<string> genelist;
        int done = 0;  
        if (g1 > 0) {            
            while (done == 0) {
                g2 = notes.find("#",g1);               
                if (g1 < 0 || g2 < 0) {
                    done = 1;
                    break;
                } else {
                    genelist.push_back(notes.substr(g1+6,g2-g1-6));                    
                    g1 = notes.find("LOCUS:",g2);
                }
            }            
        }
        if (genelist.size() > 0) {
            PROTECT(GENES = allocVector(STRSXP,genelist.size()));
            for (int s = 0;s < (int)genelist.size();s++) {   
                SET_STRING_ELT(GENES,s,mkChar(genelist.at(s).c_str()));
            }
            UNPROTECT(1);
        } else*/
        //GENES = R_NilValue;
        //SET_VECTOR_ELT(REACTION,7,GENES); SET_STRING_ELT(REACTIONNAMES,7,mkChar("genes"));

        int numOfModifiers = reactions->get(i)->getNumModifiers();
		PROTECT(GENES = allocVector(STRSXP,numOfModifiers));
		for (int m = 0;m < numOfModifiers;m++) {
			SET_STRING_ELT(GENES,m,mkChar(reactions->get(i)->getModifier(m)->getSpecies().c_str()));
		}
		SET_VECTOR_ELT(REACTION,7,GENES); SET_STRING_ELT(REACTIONNAMES,7,mkChar("genes"));



        //cout << "node level :|" << endl;
        // make the reaction edge list for convienience.
        /*SEXP EDGES, NAMES, FROM, TO, LABEL, PWAY;
        int size;
        if (LOGICAL(REVERSIBLE)[0]) size = 2*numOfReactants*numOfProducts;
        else size = numOfReactants*numOfProducts;
        PROTECT(FROM = allocVector(STRSXP,size));
        PROTECT(TO = allocVector(STRSXP,size));
        PROTECT(LABEL = allocVector(STRSXP,size));
        PROTECT(PWAY = allocVector(STRSXP,size));
        unsigned iter = 0;
        for (int r = 0;r < numOfReactants;r++) {
            for (int p = 0;p < numOfProducts;p++) {
                SET_STRING_ELT(FROM,iter,mkChar(reactions->get(i)->getReactant(r)->getSpecies().c_str()));
                SET_STRING_ELT(TO,iter,mkChar(reactions->get(i)->getProduct(p)->getSpecies().c_str()));
                SET_STRING_ELT(LABEL,iter,ST_ID);
                SET_STRING_ELT(PWAY,iter,mkChar(model->getName().c_str()));
                if (LOGICAL(REVERSIBLE)[0]) {
                        iter = iter + 1;
                        SET_STRING_ELT(FROM,iter,mkChar(reactions->get(i)->getProduct(p)->getSpecies().c_str()));
                        SET_STRING_ELT(TO,iter,mkChar(reactions->get(i)->getReactant(r)->getSpecies().c_str()));
                        SET_STRING_ELT(LABEL,iter,ST_ID);
                        SET_STRING_ELT(PWAY,iter,mkChar(model->getName().c_str()));
                }
                iter = iter + 1;
            }
        }
        PROTECT(EDGES = allocVector(VECSXP,4));
        PROTECT(NAMES = allocVector(STRSXP,4));
        SET_VECTOR_ELT(EDGES,0,FROM); SET_STRING_ELT(NAMES,0,mkChar("from"));
        SET_VECTOR_ELT(EDGES,1,TO); SET_STRING_ELT(NAMES,1,mkChar("to")); 
        SET_VECTOR_ELT(EDGES,2,LABEL); SET_STRING_ELT(NAMES,2,mkChar("label"));
        SET_VECTOR_ELT(EDGES,3,PWAY); SET_STRING_ELT(NAMES,3,mkChar("pathway"));
        setAttrib(EDGES,R_NamesSymbol,NAMES);*
        SET_VECTOR_ELT(REACTION,8,EDGES); SET_STRING_ELT(REACTIONNAMES,8,mkChar("edges"));  
		*/

        //cout << "edge level :|" << endl;
        // gets the gene association information
        
        setAttrib(REACTION,R_NamesSymbol,REACTIONNAMES);
        SET_VECTOR_ELT(REACTIONLIST,i,REACTION);

        /*SEXP cmdSexp;
        PROTECT(cmdSexp = allocVector(STRSXP, 1));
        SET_STRING_ELT(cmdSexp, 0, mkChar("Rprintf(1)");
        eval(parse(cmdSexp), R_GlobalEnv);
        //if(isNull(VECTOR_ELT(REACTIONLIST,i))){cout << "null"<<endl;}
		*/
        UNPROTECT(11);
    }
    
    setAttrib(REACTIONLIST,R_NamesSymbol,ID);
    UNPROTECT(2);
    //cout << "RacList returns :|" << endl;
    return(REACTIONLIST);
};

//
// Can be made more general if more info is present
//
SEXP getSpeciesFrame(Model *model) {
  ListOfSpecies *speciesList = model->getListOfSpecies();

  SEXP SPECIESFRAME,DIMNAMES;
  SEXP ID,NAME,COMPARTMENT,CHARGE,BOUNDARYCONDITION,CHEBI,KEGG,UNIPROT;
  
  PROTECT(ID = allocVector(STRSXP,speciesList->size()));
  PROTECT(NAME = allocVector(STRSXP,speciesList->size()));
  PROTECT(CHARGE = allocVector(INTSXP,speciesList->size()));
  PROTECT(COMPARTMENT = allocVector(STRSXP,speciesList->size()));
  PROTECT(BOUNDARYCONDITION = allocVector(LGLSXP,speciesList->size()));  
  PROTECT(CHEBI = allocVector(STRSXP,speciesList->size()));
  PROTECT(KEGG = allocVector(STRSXP,speciesList->size()));
  PROTECT(UNIPROT = allocVector(STRSXP,speciesList->size()));

  cout << "Processing All Species" << endl;
  cout << "  Number of species found: " << speciesList->size() << endl;
  for (unsigned i = 0;i < speciesList->size(); i++) {
        SET_STRING_ELT(ID,i, mkChar(speciesList->get(i)->getId().c_str()));
        SET_STRING_ELT(NAME,i, mkChar(speciesList->get(i)->getName().c_str()));
        SET_STRING_ELT(COMPARTMENT,i, mkChar(speciesList->get(i)->getCompartment().c_str()));
        INTEGER(CHARGE)[i] = speciesList->get(i)->getCharge();
        LOGICAL(BOUNDARYCONDITION)[i] = speciesList->get(i)->getBoundaryCondition();

        //Get species annotations (UniProt, ChEBI)
		int numCVTerms = speciesList->get(i)-> getNumCVTerms();
		for(int cv=0; cv<numCVTerms; cv++){
			CVTerm *cvt = speciesList->get(i)-> getCVTerm(cv);
			string URI;
			SEXP UNI_ID;
			for(int n=0;n< cvt -> getNumResources();n++){
				URI= cvt ->getResourceURI(n);
				int pos = URI.find("chebi:");
				if(pos > 0){
					PROTECT(UNI_ID = mkChar(URI.substr(pos+6).c_str()));
					SET_STRING_ELT(CHEBI,i, UNI_ID);
					UNPROTECT(1);
					continue;
				}
				pos = URI.find("kegg.compound:");
				if(pos > 0){
					PROTECT(UNI_ID = mkChar(URI.substr(pos+14).c_str()));
					SET_STRING_ELT(KEGG,i, UNI_ID);
					UNPROTECT(1);
					continue;
				}
				pos = URI.find("uniprot:");
				if(pos > 0){
					PROTECT(UNI_ID = mkChar(URI.substr(pos+8).c_str()));
					SET_STRING_ELT(UNIPROT,i, UNI_ID);
					UNPROTECT(1);
					continue;
				}
			}
  	  	}
  }  
  
  PROTECT(SPECIESFRAME = allocVector(VECSXP,8));
  PROTECT(DIMNAMES = allocVector(VECSXP,8));
  SET_VECTOR_ELT(SPECIESFRAME,0,ID);
  SET_VECTOR_ELT(SPECIESFRAME,1,NAME);
  SET_VECTOR_ELT(SPECIESFRAME,2,COMPARTMENT);
  SET_VECTOR_ELT(SPECIESFRAME,3,CHARGE);
  SET_VECTOR_ELT(SPECIESFRAME,4,BOUNDARYCONDITION);
  SET_VECTOR_ELT(SPECIESFRAME,5,CHEBI);
  SET_VECTOR_ELT(SPECIESFRAME,6,KEGG);
  SET_VECTOR_ELT(SPECIESFRAME,7,UNIPROT);

  SET_VECTOR_ELT(DIMNAMES,0,mkString("id"));
  SET_VECTOR_ELT(DIMNAMES,1,mkString("name"));
  SET_VECTOR_ELT(DIMNAMES,2,mkString("compartment"));
  SET_VECTOR_ELT(DIMNAMES,3,mkString("charge"));
  SET_VECTOR_ELT(DIMNAMES,4,mkString("boundaryCondition"));
  SET_VECTOR_ELT(DIMNAMES,5,mkString("ChEBI"));
  SET_VECTOR_ELT(DIMNAMES,6,mkString("KEGG"));
  SET_VECTOR_ELT(DIMNAMES,7,mkString("UniProt"));
  setAttrib(SPECIESFRAME,R_NamesSymbol,DIMNAMES);
  
  UNPROTECT(10);
  
  return(SPECIESFRAME);
};
