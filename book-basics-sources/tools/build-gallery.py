#!/usr/bin/env python3
"""Export the bilingual book's complete illustrations and build its figure gallery."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import html
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import pypdfium2 as pdfium
from pypdf import PdfReader

CHAPTERS = [
 ('Information and compression', 'Information et compression'),
 ('Image processing', 'Traitement d’images'),
 ('Sparsity and compressed sensing', 'Parcimonie et acquisition comprimée'),
 ('Optimal transport', 'Transport optimal'),
 ('Neural networks', 'Réseaux de neurones'),
]
TITLES = [
 ('Pixels as a message','Les pixels comme message'),
 ('Prefix codes and binary trees','Codes préfixes et arbres binaires'),
 ('Distributions and entropy','Distributions et entropie'),
 ('Coding pixel differences','Coder les différences de pixels'),
 ('The distribution of differences','La distribution des différences'),
 ('Pixels and grayscale values','Pixels et niveaux de gris'),
 ('Subsampling an image','Sous-échantillonner une image'),
 ('Quantizing grayscale values','Quantifier les niveaux de gris'),
 ('A pixel neighborhood','Le voisinage d’un pixel'),
 ('Local averaging','La moyenne locale'),
 ('Median filtering','Le filtrage médian'),
 ('Neighbors for image derivatives','Les voisins pour calculer les dérivées'),
 ('Detecting edges','Détecter les contours'),
 ('From color to grayscale','De la couleur aux niveaux de gris'),
 ('Red, green, and blue channels','Les canaux rouge, vert et bleu'),
 ('Mixing colors','Le mélange des couleurs'),
 ('Cyan, magenta, and yellow channels','Les canaux cyan, magenta et jaune'),
 ('Transforming grayscale intensity','Transformer les niveaux de gris'),
 ('Changing color contrast','Modifier le contraste en couleur'),
 ('Transposing and rotating an image','Transposer et tourner une image'),
 ('Interpolating between images','Interpoler entre deux images'),
 ('Sampling signals and images','Échantillonner les signaux et les images'),
 ('Wavelet approximation','L’approximation par ondelettes'),
 ('Blur and missing pixels','Flou et pixels manquants'),
 ('The geometry of sparsity','La géométrie de la parcimonie'),
 ('Reconstructing missing pixels','Reconstruire les pixels manquants'),
 ('Fitting data with norm constraints','Ajuster les données sous une contrainte de norme'),
 ('The single-pixel camera','La caméra à pixel unique'),
 ('Assigning bakeries to cafés','Affecter les boulangeries aux cafés'),
 ('The cost of an assignment','Le coût d’une affectation'),
 ('Matching along a metro line','L’affectation sur une ligne de métro'),
 ('Monge and planar matching','Monge et l’affectation dans le plan'),
 ('Four views of a transport plan','Quatre représentations d’un plan de transport'),
 ('Transferring a color palette','Transférer une palette de couleurs'),
 ('Barycenters of shapes','Les barycentres de formes'),
 ('Texts as word distributions','Les textes comme distributions de mots'),
 ('An image classifier','Un classifieur d’images'),
 ('Biological and artificial neurons','Neurones biologiques et artificiels'),
 ('Learning from labeled images','Apprendre à partir d’images étiquetées'),
 ('An image generator','Un générateur d’images'),
 ('Comparing point assignments','Comparer des affectations de points'),
 ('Generating and interpolating images','Générer et interpoler des images'),
]


def group(text, start):
    assert text[start] == '{'
    depth, i = 1, start + 1
    while i < len(text):
        if text[i] == '\\':
            i += 2
            continue
        if text[i] == '{': depth += 1
        elif text[i] == '}':
            depth -= 1
            if depth == 0: return text[start + 1:i], i + 1
        i += 1
    raise ValueError('Unclosed TeX group')


def labels(root, path, seen=None):
    seen = set() if seen is None else seen
    if path in seen: return {}
    seen.add(path)
    text = path.read_text()
    result = {}
    for match in re.finditer(r'\\(?:newlabel|bibcite)\{([^}]+)\}\{', text):
        body, _ = group(text, match.end()-1)
        fields, pos = [], 0
        while pos < len(body) and body[pos] == '{':
            value, pos = group(body, pos); fields.append(value)
        result[match[1]] = fields
    for child in re.findall(r'\\@input\{([^}]+)\}', text):
        result.update(labels(root, root / child, seen))
    return result


def inventory(root, lang):
    rows = []
    for path in sorted((root/'chapters-tex'/lang).rglob('*.tex')):
        text = re.sub(r'(?<!\\)%[^\n]*', '', path.read_text())
        for match in re.finditer(r'\\(?:myfigure|myfig)\s*\{|\\begin\{figure\}(?:\[[^]]*\])?', text):
            if match[0].startswith('\\begin'):
                end = text.index('\\end{figure}', match.end())
                body = text[match.end():end]
                start = body.index('\\caption{')
                caption, _ = group(body, start + len('\\caption'))
                label = re.search(r'\\label\{([^}]+)', body)[1]
                body = body[:start]
            else:
                args, pos = [], match.end()-1
                for _ in range(3):
                    value, pos = group(text, pos); args.append(value)
                    while pos < len(text) and text[pos].isspace(): pos += 1
                body, label, caption = args if match[0].startswith('\\myfig{') else [args[0], args[2], args[1]]
            rows.append({'label':label, 'body':body.strip(), 'caption':caption,
                         'source':path.relative_to(root).as_posix(),
                         'chapter':int(path.name[:2]) if path.name[:2].isdigit() else 5})
    return rows


def format_caption(caption, book_labels):
    caption = re.sub(r'\\(ref|eqref|cite)\{([^}]+)\}',
        lambda m: ('(' + book_labels[m[2]][0] + ')' if m[1]=='eqref' else
                   '[' + book_labels[m[2]][0] + ']' if m[1]=='cite' else book_labels[m[2]][0]), caption)
    caption = re.sub(r'\\norm\{([^}]+)\}', lambda m:r'\lVert '+m[1]+r'\rVert', caption)
    return ' '.join(caption.replace('~',' ').replace(r'\,','\u202f').split())


def caption_html(text):
    return ''.join('<span class="figure-math">\\('+html.escape(p[1:-1])+'\\)</span>'
        if p.startswith('$') else html.escape(p) for p in re.split(r'(\$[^$]+\$)', text))


def export(root, dest, work, item):
    folder = work/item['lang']/item['slug']
    folder.mkdir(parents=True, exist_ok=True)
    wrapper = r'\documentclass[10pt,openany]{book}'+'\n'
    if item['lang']=='fr': wrapper += r'\newcommand{\FrenchBook}{}'+'\n'
    wrapper += r'''\input{chapters-tex/shared/book-preamble}
\usepackage[active,tightpage]{preview}
\setlength{\PreviewBorder}{4pt}
\pdfinfoomitdate=1
\pdftrailerid{}
\begin{document}
\begin{preview}
\begin{minipage}{\textwidth}\centering
'''+item['body']+r'''
\end{minipage}
\end{preview}
\end{document}
'''
    source = folder/'figure.tex'; source.write_text(wrapper)
    env = dict(os.environ, TEXINPUTS=str(root)+'//:'+os.environ.get('TEXINPUTS',''))
    run = subprocess.run([shutil.which('pdflatex') or '/Library/TeX/texbin/pdflatex',
        '-halt-on-error','-interaction=nonstopmode','-recorder',f'-output-directory={folder}',str(source)],
        cwd=root, env=env, capture_output=True, text=True)
    (folder/'build.log').write_text(run.stdout)
    if run.returncode:
        raise RuntimeError(item['id']+'\n'+run.stdout[-3000:])
    log=(folder/'figure.log').read_text(errors='replace')
    problems=re.findall(r'^(?:Overfull .*|Underfull .*|Missing character:.*|.*Warning:.*)$',log,re.M)
    if problems: raise RuntimeError(item['id']+'\n'+'\n'.join(problems))
    output=dest/item['pdf'].lstrip('/');output.parent.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(folder/'figure.pdf',output)
    item['pdf_sha256']=hashlib.sha256(output.read_bytes()).hexdigest()
    print('Exported',item['lang'],item['number'],flush=True)



def share_identical_editions(items, dest):
    groups = {}
    for item in items:
        groups.setdefault((item['label'], item['jpg_sha256']), []).append(item)
    for rows in groups.values():
        if len(rows) < 2:
            continue
        first = rows[0]
        for key, extension in [('pdf', '.pdf'), ('jpg', '.jpg'), ('thumbnail', '-thumb.webp')]:
            shared = '/book-basics-sources/gallery/shared/' + first['slug'] + extension
            target = dest / shared.lstrip('/')
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(dest / first[key].lstrip('/'), target)
            for row in rows:
                (dest / row[key].lstrip('/')).unlink()
                row[key] = shared
        for row in rows:
            row['pdf_sha256'] = first['pdf_sha256']


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--book-root',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--site-root',type=Path)
    parser.add_argument('--work',type=Path)
    args=parser.parse_args();root=args.book_root.resolve()
    dest=args.site_root.resolve() if args.site_root else root.parent
    work=args.work.resolve() if args.work else root/'build/figure-gallery'
    items=[];books={};input_hashes={}
    for lang,stem in [('en','IntroImaging'),('fr','IntroImaging-FR')]:
        aux=root/'build/latex'/stem/(stem+'.aux');refs=labels(aux.parent,aux)
        pdf=root/(stem+'.pdf');reader=PdfReader(pdf)
        books[lang]={'pdf':'/book-basics-sources/'+pdf.name,'sha256':hashlib.sha256(pdf.read_bytes()).hexdigest()}
        rows=inventory(root,lang)
        rows.sort(key=lambda item:tuple(map(int,refs[item['label']][0].split('.'))))
        anchors={key for key in reader.named_destinations if key.startswith('figure.')}
        if len(rows)!=len(TITLES) or {refs[r['label']][3] for r in rows}!=anchors:
            raise ValueError('The gallery must match every numbered figure in '+lang)
        for index,row in enumerate(rows):
            number,book_page,_,anchor=refs[row['label']][:4]
            pdf_page=reader.get_destination_page_number(reader.named_destinations[anchor])+1
            if reader.page_labels[pdf_page-1]!=book_page: raise ValueError('Page label mismatch')
            slug=re.sub(r'^fig[-:]','',row['label']).replace(':','-')
            prefix=f'/book-basics-sources/gallery/{lang}/'+slug
            text=format_caption(row['caption'],refs)
            row.update(lang=lang,number=number,book_page=book_page,pdf_page=pdf_page,
                id=f'basics-{lang}-{slug}',slug=slug,title=TITLES[index][lang=='fr'],
                caption=text,caption_html=caption_html(text),pdf=prefix+'.pdf',jpg=prefix+'.jpg',
                thumbnail=prefix+'-thumb.webp',book=books[lang]['pdf'],
                chapter_title=CHAPTERS[row['chapter']-1][lang=='fr'])
            items.append(row)
    for index in range(len(TITLES)):
        en,fr=items[index],items[index+len(TITLES)]
        if en['label']!=fr['label'] or en['number']!=fr['number']: raise ValueError('Edition mismatch')
        search=' '.join([en['title'],fr['title'],en['caption'],fr['caption'],en['chapter_title'],fr['chapter_title']])
        en['search']=fr['search']=search
    with ThreadPoolExecutor(max_workers=3) as pool:
        list(pool.map(lambda item:export(root,dest,work,item),items))
    for item in items:
        with pdfium.PdfDocument(str(dest/item['pdf'].lstrip('/'))) as doc:
            if len(doc)!=1:raise ValueError('Multi-page illustration: '+item['id'])
            page=doc[0];bitmap=page.render(scale=min(2600/page.get_width(),2200/page.get_height()))
            im=bitmap.to_pil().convert('RGB');bitmap.close();page.close()
            item['width'],item['height']=im.size
            im.save(dest/item['jpg'].lstrip('/'),'JPEG',quality=95,subsampling=0,optimize=True,dpi=(300,300))
            im.thumbnail((760,560));im.save(dest/item['thumbnail'].lstrip('/'),'WEBP',quality=88,method=6)
        del item['body']
        item['jpg_sha256']=hashlib.sha256((dest/item['jpg'].lstrip('/')).read_bytes()).hexdigest()
    share_identical_editions(items, dest)
    # Input hashes make stale exports visible when an illustration or chapter changes.
    paths=set((root/'chapters-tex').rglob('*.tex'))
    for fls in work.rglob('figure.fls'):
        for name in re.findall(r'^INPUT (.+)$',fls.read_text(),re.M):
            path=Path(name);path=(root/path).resolve() if not path.is_absolute() else path.resolve()
            if root/'figures' in path.parents:paths.add(path)
    for path in paths:input_hashes[path.relative_to(root).as_posix()]=hashlib.sha256(path.read_bytes()).hexdigest()
    data={'books':books,'inputs':input_hashes,'count':len(TITLES),
          'chapters':[{'number':i+1,'title':a,'title_fr':b} for i,(a,b) in enumerate(CHAPTERS)],'figures':items}
    target=dest/'_data/basics_figures.json';target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n')
    print('Published metadata for',len(TITLES),'figures in two editions.')

if __name__=='__main__':main()
