#include <Rcpp.h>
#include <cmath>
#include <vector>

/* C++ | R INTERFACE; pruning algorithm of BM likelihood (based on diversitree:::make.bm.direct */
RcppExport SEXP bm_direct (SEXP dat, SEXP pars) 
{
	/* 
	 * all objects ordered from 1:(Nnode(phy)+Ntip(phy)) unless noted otherwise
	 * dat: a list of elements 
		* len: edge lengths (vector)
		* root: root ID
		* y: tip data 
		* order: order of internal IDs for pruning algorithm
	 * pars: rates associated with each branch 
	 */
		
	try {
		/* call in parameters associated with 'dat' object */
		Rcpp::List cache(dat);
		
		int root = Rcpp::as<int>(cache["root"]);
		int n = Rcpp::as<int>(cache["n"]);
        
        double drift = Rcpp::as<double>(cache["drift"]);

		std::vector<double> len = Rcpp::as<std::vector<double> >(cache["len"]);
		std::vector<double> y = Rcpp::as<std::vector<double> >(cache["y"]);
		std::vector<double> var = Rcpp::as<std::vector<double> >(cache["var"]);
        std::vector<double> known = Rcpp::as<std::vector<double> >(cache["given"]); //
		std::vector<int> intorder = Rcpp::as<std::vector<int> >(cache["intorder"]);
		std::vector<int> tiporder = Rcpp::as<std::vector<int> >(cache["tiporder"]);
		std::vector<int> descR = Rcpp::as<std::vector<int> >(cache["descRight"]);
		std::vector<int> descL = Rcpp::as<std::vector<int> >(cache["descLeft"]);

		std::vector<double> rates = Rcpp::as<std::vector<double> > (pars);
		
		std::vector<double> lq;
		lq.assign(n,0.0);

		double yi, ri, li, m1, m2, v1, v2, v12, m12, n12, nm, nv, m, mm, v, k;
		
		double const PIx = 4.0*atan(1.0);
		
		std::vector<double> branchinitM;
		std::vector<double> branchinitV;
		std::vector<double> branchbaseM;
		std::vector<double> branchbaseV;

		branchinitM.assign(n,0.0);
		branchinitV.assign(n,0.0);
		branchbaseM.assign(n,0.0);
		branchbaseV.assign(n,0.0);
		
		int i, z, cur, d1, d2;
		
		/* mean and variance for leaves */
		z=tiporder.size();
		for(i=0; i<z; i++){ 
			cur=tiporder[i]-1;
			yi=y[cur];
			li=len[cur];
			ri=rates[cur];
			
			branchinitM[cur] = yi;
			branchbaseM[cur] = yi + drift*li;
			branchbaseV[cur] = var[cur] + li*ri;
		}
		
		/* mean, variance, and density for edges */
		z=intorder.size();
		for(i=0; i<z; i++){ 
			cur=intorder[i]-1;
			d1=descR[cur]-1;
			d2=descL[cur]-1;
            m1=branchbaseM[d1];
			m2=branchbaseM[d2];
			v1=branchbaseV[d1];
			v2=branchbaseV[d2];
            
            v12=v1+v2;
            
            m = (((m1*v2) + (m2*v1))/v12);                  // phylogenetic mean expectation
            branchinitM[cur] = m;
            
            v = ((v1*v2)/v12);
            branchinitV[cur] = v;
            
            m12=pow((m1-m2),2);
            lq[cur] = ((-m12/(2*v12)) - (log(2*PIx*v12)/2));
            
            k=known[cur];                                   
            
            if( k == (signed)1 )                            // resolve whether node state is given (k==1)
            {
                nm=y[cur];
                nv=var[cur];
                mm=m;
                
                v12=v+nv;
                m = ((mm*nv) + (nm*v))/v12;
                branchinitM[cur] = m;

                v = (v*nv)/v12;
                branchinitV[cur] = v;
                
                m12=pow((mm-nm),2);
                lq[cur]+=((-m12/(2*v12)) - (log(2*PIx*v12)/2));
            }
            
            li=len[cur];            
            branchbaseM[cur] = m + drift*li;
            branchbaseV[cur] = v + rates[cur]*li;
		}
		
		/* compute root */ 
		cur=root-1;
		d1=descR[cur]-1;
		d2=descL[cur]-1;
		m1=branchbaseM[d1];
		m2=branchbaseM[d2];
		v1=branchbaseV[d1];
		v2=branchbaseV[d2];
		v12=v1+v2;
        
        m=(((m1*v2) + (m2*v1))/v12); //
		branchinitM[cur] = m; //
        v=((v1*v2)/v12); //
		branchinitV[cur] = v;
		m12=pow((m1-m2),2);
		lq[cur] = ((-m12/(2*v12)) - (log(2*PIx*v12)/2));
        
        /* compute root lnL (either ML or given) */
        k=known[cur];
        if(k == (signed)1 )                                
        { // given state
            nm=y[cur];                                     
            nv=var[cur] + v;
            n12=pow((nm-m),2);                             
            lq[cur] += ((-n12/(2*nv)) - (log(2*PIx*nv)/2));
        }
        else
        { // ML state
            lq[cur] += (- (log(2*PIx*v)/2));
            
        }

		
		/* PREPARE OUTPUT FOR R */
		return Rcpp::List::create(
								  Rcpp::Named("initM",branchinitM),
								  Rcpp::Named("initV",branchinitV),
								  Rcpp::Named("baseM",branchbaseM),
								  Rcpp::Named("baseV",branchbaseV),
								  Rcpp::Named("lq",lq)
								  );


    } catch( std::exception &ex ) {		
		forward_exception_to_r( ex );
    } catch(...) { 
		::Rf_error( "C++ exception: unknown reason" ); 
    }
    return R_NilValue; 
}


/* C++ | R INTERFACE; determine which branches subtended by node are not in list of excluded nodes (and their descendants)  */
RcppExport SEXP open_subtree (SEXP dat, SEXP desc) 
{
	/* 
	 * dat: a list of elements 
	 *  node: of interest
	 *  exclude: vector of excluded nodes
	 *  N: tips 
	 * desc: a list of ALL descendants from 1:(Ntip(phy)+Nnode(phy))
	 */
	
	try {
		/* call in parameters associated with 'dat' object */
		Rcpp::List cache(dat);
		Rcpp::List adesc(desc);
		int N = Rcpp::as<int>(cache["N"]);
		int n = Rcpp::as<int>(cache["n"]);
		int nn = N+n;
		int node = Rcpp::as<int>(cache["node"]);
		std::vector<int> exclude = Rcpp::as<std::vector<int> >(cache["exclude"]);
		std::vector<int> subtended;
		std::vector<int> drop;

		int i, j, k, se, sn, sd, nd, cur;

		std::vector<int> res;

		if(node>N)
		{
			/* fix exclude vector (drop) */
			se=exclude.size();
			for(i=0; i<se; i++){
				cur=exclude[i];
				drop.push_back(cur);
				if( (cur>node) & (cur<=nn) )
				{
					std::vector<int> dd=adesc[cur-1];
					sd=dd.size();
					for(j=0; j<sd; j++){
						drop.push_back(dd[j]);
					}
				}
			}
			
			std::vector<int> nodedesc = adesc[node-1];
			se=drop.size();
			if(se>0)
			{
				/* find descendants not in drop list */
				sn=nodedesc.size();
				sd=drop.size();
				for(i=0; i<sn; i++){
					cur=nodedesc[i];
					k=0;
					for(j=0; j<sd; j++){
						nd=drop[j];
						if(cur!=nd)
						{
							k+=1;
						}
					}
					if(k==sd)
					{
						res.push_back(cur);
					}
				}
			} 
			else 
			{
				sn=nodedesc.size();
				for(i=0; i<sn; i++){
					res.push_back(nodedesc[i]);
				}
			}
		}
		
		/* PREPARE OUTPUT FOR R */
		return Rcpp::wrap(res);
		
    } catch( std::exception &ex ) {		
		forward_exception_to_r( ex );
    } catch(...) { 
		::Rf_error( "C++ exception: unknown reason" ); 
    }
    return R_NilValue; 
}


RcppExport SEXP cache_descendants (SEXP phy)
{
    /* requires preorder (pruningwise ordering) of 'phylo' object */
  try {
  	Rcpp::List phylo(phy);
		int N = Rcpp::as<int>(phylo["N"]);
		int maxnode = Rcpp::as<int>(phylo["MAXNODE"]);
		std::vector<int> anc = Rcpp::as<std::vector<int> >(phylo["ANC"]);
		std::vector<int> des = Rcpp::as<std::vector<int> >(phylo["DES"]);
		
		int rows = maxnode-1;
		int root = N+1;
		
		std::vector< std::vector<int> > TIPS;
		std::vector< std::vector<int> > FDESC;
		std::vector< std::vector<int> > ADESC;
		std::vector< std::vector<int> > AANC;
        
		std::vector<int> empty;
		
		
		int i, j, k, s, t, z, dn, fd;
		
		/* initialize TIPS with known descendants (tips and root), otherwise leave empty */
		std::vector<int> cur;
		for(i = 0; i < maxnode; i++) {
			FDESC.push_back(empty);
            AANC.push_back(empty);
			if(i < N)
			{
				cur.push_back(i+1);
				TIPS.push_back(cur);
				ADESC.push_back(cur);
			}
			else
			{
				TIPS.push_back(empty);
				ADESC.push_back(empty);
			}
			cur.clear();
		}
		
		/* store nodes associated with root -- TIPS */
		for(i=0; i<N; i++){
			cur.push_back(i+1);
		}
		TIPS.at(N)=cur;
		
		/* store nodes associated with root -- ALL */
		cur.clear();
		for(i=0; i<maxnode; i++){
			if(i!=N)
			{
				cur.push_back(i+1);
			}
		}
		ADESC.at(N)=cur;
		
		/* store nodes associated with root -- FIRST */
		cur.clear();
		for(i=0; i<rows; i++){
			int idx = anc.at(i);
			if(idx==root)
			{
				cur.push_back(des.at(i));
			}
		}
		FDESC.at(N)=cur;
		
		/* collect descendants of each node in edge matrix (using pruningwise order to eliminate unnecessary computations) */
		for(i = 0; i < rows; i++){
			int nd = des.at(i);
			if(nd > N) {
				/* find descendant nodes */
				std::vector<int> subtends;
				for(j=0; j<rows; j++){
					int idx = anc.at(j);
					if(idx==nd)
					{
						subtends.push_back(des.at(j));
					}
				}
				FDESC.at(nd-1)=subtends;
				
				/* find immediate descendants of nd */
				std::vector<int> subtendedtips;
				std::vector<int> subtendednodes;
				s=subtends.size();
				for(k = 0; k < s; k++){
					
					/* find nodes subtended by immediate descendants of nd */
					fd = subtends.at(k);
					subtendednodes.push_back(fd);
                    
					if(fd<root) {
						subtendedtips.push_back(fd);
					} else {
						std::vector<int> descnodes = ADESC[fd-1];
						t=descnodes.size();
						for(z = 0; z < t; z++){
							dn=descnodes[z];
							if(dn<root)
							{
								subtendedtips.push_back(dn);
							}
							subtendednodes.push_back(dn);
						}
					}
				}
				/* store tips associated with nd into main list */
				TIPS.at(nd-1)=subtendedtips;
				ADESC.at(nd-1)=subtendednodes;
                s=subtendednodes.size();
                
                /* store ancestors */
                for(k=0; k<s; k++){
                    int idx = subtendednodes.at(k);
                    std::vector<int> ancnodes = AANC[idx-1];
                    ancnodes.push_back(nd);
                    AANC.at(idx-1)=ancnodes;
                }
			}
		}
		
		for(i=0; i<N; i++){
			ADESC.at(i)=empty;
		}
        for(k=0; k<maxnode; k++){
            int idx = k+1;
            if(idx!=root){
                std::vector<int> ancnodes = AANC[k];
                ancnodes.push_back(root);
                AANC.at(k)=ancnodes;
            }
        }
		return Rcpp::List::create(Rcpp::Named("tips",TIPS),
								  Rcpp::Named("fdesc",FDESC),
								  Rcpp::Named("adesc",ADESC),
                                  Rcpp::Named("anc",AANC));
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "C++ exception: unknown reason" );
	}
	return R_NilValue;
}
