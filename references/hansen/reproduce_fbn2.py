"""Reproducible Fbn2 relative-MSD analysis. Uses NumPy and h5py, no fitted model.

Run: python reproduce_fbn2.py --data /path/to/downloads --output /path/to/derived
Input files are checksum-verified Zenodo downloads listed in download manifests.
"""
from pathlib import Path
import csv, hashlib, json, collections, argparse
import numpy as np
import h5py

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--data',type=Path,required=True,help='Directory containing the downloaded source files; hashes are checked against the manifest beside this script')
parser.add_argument('--output',type=Path,required=True,help='Directory for derived curves and catalog')
args=parser.parse_args()
DATA=args.data.resolve()
OUT=args.output.resolve()
OUT.mkdir(parents=True,exist_ok=True);(OUT/'curves').mkdir(exist_ok=True)
FRAME_S=20; MAX_LAG=180; BOOT=1000; SEED=20261102
NOISE_SIGMA_UM=.203
NOISE_OFFSET=2*NOISE_SIGMA_UM**2
BOOT_URL='https://github.com/ahansenlab/polysim_bootcamp/blob/7f273ab648e7c606c7212702bb5caa1604810006/tutorial/calibration.ipynb'
SOURCE='https://zenodo.org/records/5770531'
PAPER='https://doi.org/10.1126/science.abn6583'
MANIFEST=json.loads((Path(__file__).resolve().parent/'download-manifest.json').read_text())
MANIFEST={x['file']:x for x in MANIFEST}
for filename in list(k for k in MANIFEST if k.endswith('.tagged_set.tsv'))+['2022_Fbn2_full_data.h5']:
    if hashlib.sha256((DATA/filename).read_bytes()).hexdigest()!=MANIFEST[filename]['sha256']:
        raise ValueError(f'Source checksum mismatch: {filename}')
LABELS={'C27':'C27 · ΔTAD (10 kb label separation)','C36':'C36 · wild type','C65':'C65 · deleted CTCF sites','CTCF_0_hr':'CTCF degron · untreated','CTCF_2_hr':'CTCF depletion · 2 h','CTCF_4_hr':'CTCF depletion · 4 h','Rad21_0_hr':'RAD21 degron · untreated','Rad21_2_hr':'RAD21 depletion · 2 h','Rad21_4_hr':'RAD21 depletion · 4 h','WAPL_0_hr':'WAPL degron · untreated','WAPL_4_hr':'WAPL depletion · 4 h','WAPL_6_hr':'WAPL depletion · 6 h'}
SLUGS={'C27':'c27','C36':'c36','C65':'c65',**{f'{p}_{h}_hr':f'{p.lower()}-{h}h' for p,hs in [('CTCF',[0,2,4]),('Rad21',[0,2,4]),('WAPL',[0,4,6])] for h in hs}}


def dump(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')

def load_data():
    groups={};checks={}
    with h5py.File(DATA/'2022_Fbn2_full_data.h5') as f:
        for i in sorted(f['data/_data'],key=int):
            tags=[v for k,v in f['data/_tags/'+i].attrs.items() if not k.startswith('_')]
            cond=next(t for t in tags if t in LABELS)
            g=f['data/_data/'+i];r=np.asarray(g['data'][0],float)
            assert r.ndim==2 and r.shape[1]==3
            groups.setdefault(cond,[]).append({'relative':r,'id':i,'movie':None})
    for file in sorted(DATA.glob('*.tagged_set.tsv')):
        cond=file.name.split('.')[0]; records=collections.defaultdict(list)
        for row in csv.DictReader(file.open(),delimiter='\t'):records[int(row['id'])].append(row)
        actual=[];maxerr=0.;dist_err=0; duplicate_rows=0
        assert len(records)==len(groups[cond])
        for rows, h in zip([records[k] for k in sorted(records)],groups[cond]):
            unique={int(row['t']):row for row in rows}; duplicate_rows+=len(rows)-len(unique); rows=list(unique.values())
            frames=np.array([int(x['t']) for x in rows]);assert len(np.unique(frames))==len(frames)
            xyz=np.array([[float(x[k]) for k in ['x','y','z','x2','y2','z2']] for x in rows]);r=xyz[:,3:]-xyz[:,:3]
            assert np.isfinite(r).all()
            vals=h['relative'][frames-frames.min()]
            assert np.isfinite(h['relative']).all(1).sum()==len(frames)
            err=float(np.max(abs(vals-r)));maxerr=max(maxerr,err);assert err<2e-7
            dist=np.array([float(x['dist']) for x in rows]);dist_err=max(dist_err,float(np.max(abs(dist-np.linalg.norm(r,axis=1)))))
            movies={x['movie_index'] for x in rows};assert len(movies)==1
            full=np.full((frames.max()-frames.min()+1,3),np.nan);full[frames-frames.min()]=r
            actual.append({'relative':full,'id':rows[0]['id'],'movie':movies.pop()})
        groups[cond]=actual;checks[cond]={'h5_tsv_max_coordinate_difference_um':maxerr,'provided_distance_max_difference_um':dist_err,'tracks_match':len(actual),'duplicate_frames_resolved_last_row_to_match_author_h5':duplicate_rows}
    return groups,checks


def summarize(tracks,cond,split,noise=False):
    # WAPL6h lacks acquisition labels; tracks remain separate, never falsely called movies.
    has_movie=all(x['movie'] is not None for x in tracks)
    clusters=sorted({x['movie'] if has_movie else x['id'] for x in tracks})
    idx={c:i for i,c in enumerate(clusters)}; n=len(clusters)
    sums=np.zeros((n,MAX_LAG));counts=np.zeros((n,MAX_LAG),dtype=np.int64)
    dist_sums=np.zeros((n,3));dist_n=np.zeros(n,dtype=np.int64);raw_dist=[]
    durations=[];nvalid=[]
    for tr in tracks:
        a=tr['relative'];i=idx[tr['movie'] if has_movie else tr['id']];v=np.isfinite(a).all(1);d=np.linalg.norm(a[v],axis=1)
        raw_dist.extend(d.tolist());dist_sums[i]+=[d.sum(),np.square(d).sum(),np.sum(d<.2)];dist_n[i]+=len(d)
        durations.append((len(a)-1)*FRAME_S);nvalid.append(int(v.sum()))
        for lag in range(1,min(MAX_LAG,len(a)-1)+1):
            good=v[lag:]&v[:-lag];x=a[lag:][good]-a[:-lag][good]
            sums[i,lag-1]+=np.sum(x*x);counts[i,lag-1]+=len(x)
    total_count=counts.sum(0);assert (total_count>0).all()
    raw_msd=sums.sum(0)/total_count;offset=NOISE_OFFSET if noise else 0.;msd=raw_msd-offset
    assert (msd>0).all(),(cond,split,msd.min())
    # Available acquisition (movie_index), not a claim of independent biological replicate.
    rng=np.random.default_rng(SEED)
    weights=rng.multinomial(n,np.full(n,1/n),size=BOOT)
    boot_count=weights@counts;boot_sum=weights@sums
    bmsd=np.divide(boot_sum,boot_count,out=np.full_like(boot_sum,np.nan),where=boot_count>0)-offset
    lo,hi=np.nanquantile(bmsd,[.025,.975],axis=0)
    bdist=(weights@dist_sums)/(weights@dist_n)[:,None];dlo,dhi=np.quantile(bdist,[.025,.975],axis=0)
    dm=dist_sums.sum(0)/dist_n.sum();lags=np.arange(1,MAX_LAG+1)*FRAME_S
    window=(lags>=60)&(lags<=1200)
    slope=float(np.polyfit(np.log(lags[window]),np.log(msd[window]),1)[0])
    hist,edges=np.histogram(raw_dist,bins=np.linspace(0,3,61)); overflow=int(np.sum(np.array(raw_dist)>3))
    warnings=['MSD displacement pairs overlap; counts are not numbers of independent observations.',
              'Movie-index clusters can contain shared biological preparation; intervals quantify acquisition-level resampling only.',
              'No biological clock has been assigned to a Chromatin Studio model by this reference.',
              'The pair-vector observable cancels common translation but does not remove nuclear rotation.']
    if noise:warnings+=['Noise subtraction assumes independent frame-localization error with the tutorial fixed parameter; uncertainty in that parameter and exposure blur are not included.']
    else:warnings+=['Localization noise remains in this curve. Compare an observation-matched simulation or supply a condition-specific noise model before clock fitting.']
    if not has_movie:warnings+=['WAPL 6h is available only in the H5, which omits movie IDs: uncertainty is track-level and no acquisition-held-out split is possible.']
    source_file=cond+'.tagged_set.tsv' if has_movie else '2022_Fbn2_full_data.h5'
    manifest=MANIFEST[source_file]
    analysis={'dataset':'Gabriele et al. 2022 Fbn2 paired-locus tracking','condition':cond,'condition_title':LABELS[cond],
      'species':'Mus musculus','cell_type':'mouse embryonic stem cells','cell_state':'G1/early S-phase selected in source study',
      'locus':'Fbn2','genomic_pair_separation_bp':10000 if cond=='C27' else 515000,
      'source_publication':PAPER,'source_dataset':SOURCE if has_movie else 'https://zenodo.org/records/14600051',
      'source_file':source_file,'source_sha256':manifest['sha256'],'source_license':'CC-BY-4.0',
      'processing':'QC-filtered coordinates; duplicate frames use the last source row, verified against author H5;  R(t)=r_mScarlet(t)-r_EGFP(t); average |R(t+lag)-R(t)|^2 over all valid exact-frame-lag pairs; no interpolation across gaps.',
      'estimator':'pooled valid displacement-pair weighting','half_relative_msd':False,'frame_interval_s':FRAME_S,
      'split':split,'split_unit':'movie_index' if has_movie else 'unavailable',
      'cluster_ids':clusters,'tracks':len(tracks),'clusters':n,'coordinate_observations':int(dist_n.sum()),
      'track_duration_s':{'minimum':min(durations),'median':float(np.median(durations)),'maximum':max(durations)},
      'valid_frames_per_track':{'minimum':min(nvalid),'median':float(np.median(nvalid)),'maximum':max(nvalid)},
      'maximum_lag_s':MAX_LAG*FRAME_S,'maximum_lag_rule':'At most half the nominal 365-frame recording duration; pairs are included only where both measured endpoints exist.',
      'noise':{'method':'fixed_additive_MSD_offset' if noise else 'none','relative_vector_error_parameter_um':NOISE_SIGMA_UM if noise else None,
               'subtracted_msd_um2':offset,'parameter_source':BOOT_URL if noise else None,
               'condition_scope':'RAD21 2h only' if noise else None,'raw_msd_um2':raw_msd.tolist() if noise else None},
      'uncertainty':{'method':'percentile cluster bootstrap','cluster_unit':'movie_index' if has_movie else 'trajectory',
        'resamples':BOOT,'seed':SEED,'confidence':.95,'lower_msd_um2':lo.tolist(),'upper_msd_um2':hi.tolist(),
        'parameter_uncertainty_included':False,'biological_replicate_independence_established':False},
      'distance_statistics':{'coordinate_observations':int(dist_n.sum()),'mean_distance_um':float(dm[0]),'mean_distance_95pct_interval_um':[float(dlo[0]),float(dhi[0])],
        'mean_squared_distance_um2':float(dm[1]),'mean_squared_distance_95pct_interval_um2':[float(dlo[1]),float(dhi[1])],
        'twice_mean_squared_distance_um2':float(2*dm[1]),'twice_mean_squared_distance_95pct_interval_um2':[float(2*dlo[1]),float(2*dhi[1])],
        'noise_corrected_mean_squared_distance_um2':float(dm[1]-offset/2) if noise else None,
        'noise_corrected_twice_mean_squared_distance_um2':float(2*dm[1]-offset) if noise else None,
        'noise_corrected_twice_mean_squared_distance_95pct_interval_um2':[float(2*dlo[1]-offset),float(2*dhi[1]-offset)] if noise else None,
        'moment_caveat':'2 E[|R|^2] is the tutorial second-moment proxy, not a directly observed long-lag plateau; equality requires stationary decorrelated zero-mean pair vectors.',
        'histogram':{'edges_um':edges.tolist(),'counts':hist.tolist(),'above_last_edge':overflow}},
      'suggested_fit_window':{'minimum_s':60,'maximum_s':1200,'points':int(window.sum()),'minimum_pairs_per_lag':int(total_count[window].min()),
        'maximum_pairs_per_lag':int(total_count[window].max()),'descriptive_log_log_slope':slope,
        'selection':'Fixed exploratory window, excludes first two lags and long-lag tail; model stationarity, sensitivity and shape agreement must still be checked.'}}
    return {'lag':lags.tolist(),'msd':msd.tolist(),'counts':total_count.tolist(),'dimension':3,'observable':'relative_pair',
            'drift_correction':'none','lag_unit':'s','msd_unit':'um2','source':f'Gabriele et al. 2022; Fbn2 {LABELS[cond]}; {split}; '+('tutorial noise correction; ' if noise else 'uncorrected localization noise; ')+analysis['source_dataset'],
            'analysis':analysis,'warnings':warnings}


def main():
    groups,checks=load_data();references=[];condition_summaries=[]
    for cond,tracks in sorted(groups.items()):
        slug=SLUGS[cond];has_movie=all(t['movie'] for t in tracks)
        # Split baseline by acquisition before calculating any metric. Validation conditions remain held out in entirety.
        if cond=='Rad21_2_hr':
            movies=sorted({t['movie'] for t in tracks},key=lambda m:hashlib.sha256(('Fbn2-v1:'+m).encode()).hexdigest())
            holdout=set(movies[:max(1,int(round(.2*len(movies))))]);parts={'fit':[t for t in tracks if t['movie'] not in holdout],
                'holdout':[t for t in tracks if t['movie'] in holdout],'all':tracks}
        else:parts={'validation':tracks}
        for split,selected in parts.items():
            for noise in ([False,True] if cond=='Rad21_2_hr' else [False]):
                variant='corrected' if noise else 'raw';rid=f'fbn2-{slug}-{split}-{variant}'
                curve=summarize(selected,cond,split,noise);file=f'curves/{rid}.json';dump(OUT/file,curve)
                is_default=(cond=='Rad21_2_hr' and split=='fit' and noise) or cond!='Rad21_2_hr'
                ref={'id':rid,'title':f'Fbn2 · {LABELS[cond]} · {split} · '+('noise corrected' if noise else 'raw MSD'),
                     'condition':cond,'condition_title':LABELS[cond],'role':'baseline_fit' if split=='fit' else ('baseline_holdout' if split=='holdout' else ('reference_only' if split=='all' else 'perturbation_validation')),
                     'is_default_for_condition':is_default,'recommended':cond=='Rad21_2_hr' and split=='fit' and noise,
                     'curve_file':file,'source_kind':'experimental','dimension':3,'observable':'relative_pair','drift_correction':'none',
                     'default_fit_min_s':60,'default_fit_max_s':1200,'genomic_pair_separation_bp':curve['analysis']['genomic_pair_separation_bp'],
                     'noise_corrected':noise,'tracks':len(selected),'clusters':curve['analysis']['clusters'],'points':len(curve['lag']),
                     'source_url':curve['analysis']['source_dataset'],'publication_url':PAPER,'license':'CC-BY-4.0',
                     'notes':'Noise parameter is fixed from the Hansen tutorial, not estimated from this subset.' if noise else 'Localization noise is retained; this is a structural/dynamical validation reference, not a ready-to-use biological calibration.'}
                references.append(ref)
                if is_default:condition_summaries.append({'condition':cond,'reference_id':rid,'tracks':len(tracks),'analysis':curve['analysis']})
    references.sort(key=lambda r:(not r['recommended'],not r['is_default_for_condition'],r['condition'],r['id']))
    catalog={'version':1,'title':'Hansen lab chromatin calibration references','source_scope':'Gabriele2022 all 12 conditions; 11 original filtered TSVs plus WAPL6h from the Jusuf2026 H5 re-deposit.',
             'default_reference_id':'fbn2-rad21-2h-fit-corrected','status':'empirical_references_ready_no_production_clock_fitted',
             'references':references,'warnings':['Different cell lines, perturbations and label spacings are not pooled.','These experimental references do not themselves calibrate any polymer or protein parameter.','Only RAD21 2h is designated baseline fitting data; other conditions and its acquisition holdout are reserved for validation.']}
    dump(OUT/'catalog.json',catalog);dump(OUT/'analysis-report.json',{'version':1,'created_date':'2026-11-02','conditions':condition_summaries,'h5_tsv_cross_checks':checks,
       'baseline_split':{'method':'sort movie IDs by SHA-256 of Fbn2-v1: plus movie ID; first rounded 20% held out before analysis','seed_is_implicit_in_fixed_prefix':True},
       'processing_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    dump(OUT/'download-manifest.json',list(MANIFEST.values()))
    for r in references:
        if r['is_default_for_condition']:
            c=json.loads((OUT/r['curve_file']).read_text());a=c['analysis'];print(r['id'],a['tracks'],a['clusters'],round(a['distance_statistics']['mean_distance_um'],4),a['suggested_fit_window'])

if __name__=='__main__':main()
