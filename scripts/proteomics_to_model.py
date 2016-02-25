class Prot2Model(object):
    
    def __init__(self):
        self.genes_to_active_reactions = self
   
    def filter_idle_reactions(self.V,model,use_cache=True):
        """
           This function is important for dealing with enzymes that carry several
           metabolic reactions. In such cases, we consider that the enzyme is 
           evenlly distributed between all reactions, i.e., if there are 1000 copies
           of the enzyme and it carries 4 reactions in a given condition, each reaction 
           is catalyzed by 250 copies of the enzyme. 
           
           To do so, in each condition we ask what are the flux carrying reactions
           according to the FBA solution. Then we iterate over all genes in the
           model (enzymes) and for each gene we take the intercept of the reactions
           it can carry with flux carrying reactions from FBA.
           
              
            
        """
       
        if use_cache:        
            with open('../cache/genes_to_flux_carrying_reactions.p', 'rb') as fp:
                return pickle.load(fp)
        out = {}
        for c in V.columns:
            out[c] = {}        
            v = V[c].dropna()
            for g in model.genes:
                rxns = {r.id for r in list(g.reactions)} & set(v.index)
                if len(rxns)>0:
                    out[c][g.id] = rxns

        with open('../cache/genes_to_flux_carrying_reactions.p', 'wb') as fp:
            pickle.dump(out, fp)
        return out
        
    def flux_carrying_reactions_to_enzymes(V,E,model,use_cache=True):
        if use_cache:       
            with open('../cache/flux_carrying_reactions_to_enzymes.p', 'rb') as fp:
                return pickle.load(fp)
        try:
            V = V.drop('flux_counter')
        except ValueError:
            print "flux couter already removed"
        mapper = {}
        for c in V.columns:
            mapper[c] = {}
            #use only flux carrying reactions in a given condition
            vc = V[c].dropna()
            for rid in vc.index:
                r = model.reactions.get_by_id(rid)
                genes = {g.id:g for g in r.genes}
                # annoing gene in the model - just ignore the reaction it carries
                if 's0001' in genes: continue
                mapper[c][r.id] = {}
                for i, (gid, g) in enumerate(genes.iteritems()):
                    rxns = {r.id for r in list(g.reactions)} & set(vc.index)
                    mapper[c][rid][gid] = float(len(rxns))
        with open('../cache/flux_carrying_reactions_to_enzymes.p', 'wb') as fp:
            pickle.dump(mapper, fp)
        return mapper
