#!/usr/bin/env python3
"""Inventory every displayed figure in the active book, including nested inputs."""
from pathlib import Path
import json
import re
import unicodedata
from figure_tex import braced, book_labels
from build_book import chapter_list

ROOT = Path(__file__).resolve().parents[1]

def active(text):
    text = re.sub(r'(?<!\\)%[^\n]*', lambda m: ' '*len(m[0]), text)
    return re.sub(r'\\(?:iffalse\b|if\s+0\b).*?\\fi\b', lambda m: ''.join('\n' if c=='\n' else ' ' for c in m[0]), text, flags=re.S)

def slug(text):
    text = text.replace(r'\&', ' and ')
    text = unicodedata.normalize('NFKD', text).encode('ascii','ignore').decode().lower()
    return re.sub(r'[^a-z0-9]+','-',text).strip('-')

def graphics(text):
    results=[]
    pattern = re.compile(r'\\(includegraphics|image|BookDrawing|BookPanel)\b(?:\[[^\]]*\])?')
    for m in pattern.finditer(text):
        count=3 if m[1] in ('image','BookPanel') else 1
        pos=m.end(); args=[]
        for _ in range(count):
            value,pos=braced(text,pos);args.append(value)
        name=args[0]+'/'+args[2] if m[1]=='image' else args[-1]
        base=ROOT/'figures'/name
        candidates=[base] if base.suffix else [base.with_suffix(ext) for ext in ('.pdf','.png','.jpg','.jpeg','.eps')]
        if base.suffix=='.eps':candidates=[base.with_name(base.stem+'-eps-converted-to.pdf'),base]
        found=next((p for p in candidates if p.is_file()),None)
        results.append({'reference':name,'path':str(found.relative_to(ROOT)) if found else None,
                        'command':m[1],'args':args,'start':m.start(),'end':pos})
    return results

def figures(path, chapter):
    raw=path.read_text();text=active(raw); result=[]
    pattern=re.compile(r'\\(myfigure|wrapf|wrapfSimple|input)\b|\\begin\{figure\*?\}')
    cursor=0
    while (m:=pattern.search(text,cursor)):
        kind=m[1] or 'figure'; pos=m.end()
        if kind=='input':
            name,end=braced(text,pos);p=ROOT/name
            if not p.suffix:p=p.with_suffix('.tex')
            result.extend(figures(p,chapter));cursor=end;continue
        args=[]
        if kind!='figure':
            for _ in range({'myfigure':3,'wrapf':2,'wrapfSimple':1}[kind]):
                arg,pos=braced(text,pos);args.append(arg)
            end=pos
            body=args[0] if kind=='myfigure' else r'\includegraphics[width=.5\linewidth]{'+args[0]+'}'
            caption=args[1] if len(args)>1 else ''
            label=args[2].strip() if kind=='myfigure' else None
            continued=False
        else:
            ending=re.search(r'\\end\{figure\*?\}',text[pos:]);assert ending,path
            end=pos+ending.end();body=text[pos:pos+ending.start()]
            body=re.sub(r'^\[[^\]]*\]','',body)
            cap=re.search(r'\\caption(?:\[[^\]]*\])?',body)
            if cap:
                caption,capend=braced(body,cap.end());body=body[:cap.start()]+body[capend:]
            else:caption=''
            matches=re.findall(r'\\label\{([^}]+)\}',body+'\n'+caption);label=matches[0] if matches else None
            continued=r'\BookContinuedFigure' in body
            body=re.sub(r'\\label\{[^}]+\}|\\BookContinuedFigure','',body)
        caption=re.sub(r'\\label\{[^}]+\}', '', caption)
        context_start=max(0,text.rfind('\\section',0,m.start()),m.start()-4500)
        context=text[context_start:min(len(text),end+2000)]
        result.append({'chapter':chapter,'source':str(path.relative_to(ROOT)),'start':m.start(),'end':end,
            'line':raw.count('\n',0,m.start())+1,'kind':kind,'body':body.strip(),'caption_source':caption.strip(),
            'book_label':label,'continued':continued,'assets':graphics(body),'context':context})
        cursor=end
    return result

def inventory():
    labels=book_labels(ROOT,ROOT/'build',chapter_list());result=[]
    for chno,chapter in enumerate(chapter_list(),1):
        path=ROOT/'chapters'/f'{chapter}.tex'; t=active(path.read_text())
        title,_=braced(t,t.index(r'\chapter')+len(r'\chapter'))
        found=figures(path,chapter);number=0
        aux=(ROOT/'build/book/chapters'/f'{chapter}.aux').read_text()
        caption_refs={}
        for line in aux.splitlines():
            if line.startswith(r'\@writefile{lof}{\contentsline {figure}'):
                body,_=braced(line,len(r'\@writefile{lof}'))
                pos=body.index(r'\numberline')+len(r'\numberline')
                num,pos=braced(body,pos);cap,pos=braced(body,pos)
                outer_end=pos+1
                page,outer_end=braced(body,outer_end);anchor,_=braced(body,outer_end)
                caption_refs[num]=(page,re.sub(r'\\relax\s*$', '', cap.replace(r'\ignorespaces ','')).strip(),anchor)
        for item in found:
            if item['caption_source'] and not item['continued']: number+=1
            num=f'{chno}.{number}'
            info=labels.get(item['book_label'])
            if info: assert info[0]==num,(chapter,num,item['book_label'],info)
            if item['caption_source']:
                page,caption,anchor=(info[1],info[2],info[3]) if info else caption_refs[num]
                if not info:caption=item['caption_source']
            else:
                page,caption,anchor=None,'Unnumbered illustration',None
            stable=item['book_label'] or (f'unnumbered-{len(result)}' if not item['caption_source'] else 'figure-'+num)
            stable=re.sub(r'^(?:fig[:\-])?(?:review-)?','',stable)
            if stable.startswith('figure-'):
                stable=slug(re.sub(r'\$[^$]*\$|\\[A-Za-z]+','',caption))[:70]
            item.update(chapter_title=title,chapter_slug='chapter-'+slug(title),figure_number=num if item['caption_source'] else None,
                book_page=page,book_anchor=anchor,caption_latex=caption,figure_slug=slug(stable),
                id=chapter+'--'+slug(stable))
            result.append(item)
    return result

if __name__=='__main__':
    result=inventory();out=ROOT/'build/figure-regeneration';out.mkdir(parents=True,exist_ok=True)
    (out/'inventory.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    lines=[]
    for i in result:
        lines.append(f"{i['figure_number'] or 'unnumbered'} | {i['id']} | {i['caption_source']}\n  "+', '.join(a['reference'] for a in i['assets']))
    (out/'inventory.txt').write_text('\n'.join(lines)+'\n')
    print(f'{len(result)} figure displays; {sum(len(i["assets"]) for i in result)} asset references.')
    missing=[a for i in result for a in i['assets'] if a['path'] is None]
    print('Unresolved assets:',missing)
