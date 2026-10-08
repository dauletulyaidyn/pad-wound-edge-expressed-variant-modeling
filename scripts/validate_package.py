from pathlib import Path
import ast, csv, hashlib, json
ROOT = Path(__file__).resolve().parents[1]
def main():
    manifest=json.loads((ROOT/'provenance.json').read_text())
    for item in manifest:
        p=ROOT/item['file']
        assert p.is_file(), p
        assert hashlib.sha256(p.read_bytes()).hexdigest()==item['sha256'], p
    for p in (ROOT/'scripts').rglob('*.py'):
        ast.parse(p.read_text(encoding='utf-8-sig'), filename=str(p))
    with (ROOT/'supplementary/tables/Table_S1_run_metadata_reconciliation.csv').open(newline='',encoding='utf-8') as f:
        rows=list(csv.DictReader(f))
    assert len(rows)==14 and len({r['gsm'] for r in rows})==7
    assert sum(r['source_condition']=='WE' for r in rows)==6
    assert sum(r['source_condition']=='UWE' for r in rows)==8
    conflicts=sum(r['labels_agree']=='false' for r in rows)
    assert conflicts>0
    print(json.dumps({'integrity':'PASS','source_files':len(manifest),'runs':14,'biological_samples':7,'label_conflicts':conflicts,'manuscript_results':'NOT_REPRODUCED'},indent=2))
if __name__=='__main__': main()
