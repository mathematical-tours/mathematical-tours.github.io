#!/usr/bin/env python3
"""Derive the MyST reading edition from the active LaTeX, without editing it.

Pandoc expands the book's mathematical macros and parses prose. Balanced TeX
preprocessing preserves the publication's counters, figures and theorem titles.
The reference snapshot permits a clean checkout to build without TeX Live.
"""
from pathlib import Path
import sys, re, json, subprocess, hashlib, html, shutil, collections, functools
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'tools'))
from figure_tex import braced
from figure_inventory import active, figures, graphics, slug
from build_book import chapter_list
WEB = ROOT/'myst'
THEOREMS = {'thm':'Theorem','prop':'Proposition','lem':'Lemma','cor':'Corollary','defn':'Definition','alg':'Procedure','rem':'Remark','exmp':'Example','ass':'Main assumption'}
REPORT = {'chapters':[], 'unknown_tex':[], 'references':[], 'equations':[], 'figures':[], 'warnings':[]}
LABELS = {}; OWNERS = {}; ALIASES = {}; PAGE = ''; CHNO = 0


def web_id(label):
    # Apostrophes are typographically transformed by Markdown. Stable web IDs
    # avoid that ambiguity without renaming the canonical LaTeX label.
    return label.replace("'", '-').replace('’', '-')


def references():
    snapshot = WEB/'reference-metadata.json'
    if (ROOT/'build/book/FundationsDataScience.aux').exists():
        labels={}; owners={}; bib={}
        for chapter in ['FundationsDataScience','abstract',*chapter_list()]:
            p=ROOT/'build/book'/('' if chapter=='FundationsDataScience' else 'chapters-tex')/(chapter+'.aux')
            for line in p.read_text().splitlines():
                if line.startswith(r'\newlabel{'):
                    key,pos=braced(line,len(r'\newlabel'));body,_=braced(line,pos);values=[];pos=0
                    while pos<len(body): value,pos=braced(body,pos);values.append(value)
                    labels[key]=values;owners[key]=chapter
                if line.startswith(r'\bibcite{'):
                    key,pos=braced(line,len(r'\bibcite'));value,_=braced(line,pos);bib[key]=value
        data={'labels':labels,'owners':owners,'bibliography':bib}
        snapshot.write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n')
    return json.loads(snapshot.read_text())


def macro_definitions():
    text=active((ROOT/'chapters-tex/shared/mystyle.sty').read_text())
    chunks=[]
    for m in re.finditer(r'\\(?:newcommand|renewcommand|DeclareMathOperator)\*?\s*\{',text):
        name,p=braced(text,m.end()-1)
        if name in [r'\myfigure',r'\wrapf',r'\wrapfSimple',r'\image']:continue
        arg=re.match(r'\s*\[\d+\]',text[p:])
        if arg:p+=arg.end()
        try:_,end=braced(text,p)
        except RuntimeError:continue
        chunks.append(text[m.start():end])
    # Browser math equivalents for purely typographic TeX commands.
    chunks.extend([r'\newcommand{\ones}{\mathbb{1}}',r'\newcommand{\zeros}{\mathbf{0}}',
      r'\renewcommand{\eqdef}{\mathrel{:=}}',r'\renewcommand{\Lun}{L^1}',r'\renewcommand{\Ldeux}{L^2}',r'\renewcommand{\Linf}{L^\infty}',
      r'\newcommand{\Symb}{s}',r'\newcommand{\noise}{w}',r'\newcommand{\piN}{\frac{2\mathrm{i}\pi}{N}}',
      r'\newcommand{\bk}[1]{\left[#1\right]}',r'\newcommand{\pdb}[2]{\bk{\pd{#1}{#2}}}'])
    return '\n'.join(chunks)+'\n'
PREAMBLE=macro_definitions()


@functools.lru_cache(maxsize=4096)
def pandoc(source):
    # These chapter-local definitions are also needed when rendering captions
    # separately. Keep their local override instead of declaring them twice.
    source=re.sub(r'\\newcommand(?=\s*\{\\(?:piN|Symb|noise|bk|pdb)\})',r'\\renewcommand',source)
    result=subprocess.run(['pandoc','-f','latex+raw_tex','-t','json'],input=PREAMBLE+source,text=True,capture_output=True,cwd=ROOT)
    if result.returncode:raise RuntimeError(result.stderr)
    if result.stderr: REPORT['warnings'].append(result.stderr)
    return json.loads(result.stdout)['blocks']


def replace_commands(text,name,callback,nargs=1):
    out=[];start=0
    for m in re.finditer(r'\\'+name+r'\b',text):
        if m.start()<start:continue
        p=m.end();args=[]
        for _ in range(nargs):a,p=braced(text,p);args.append(a)
        out.extend([text[start:m.start()],callback(*args)]);start=p
    out.append(text[start:]);return ''.join(out)


def split_rows(body):
    out=[];depth=0;start=0;i=0;env=0
    while i<len(body):
        if body[i]=='\\':
            if body[i:i+2]=='\\\\' and depth==0 and env==0:
                out.append(body[start:i]);i+=2
                gap=re.match(r'\s*\[[^\]]*\]',body[i:])
                if gap:i+=gap.end()
                start=i;continue
            nested=re.match(r'\\(begin|end)\{[^}]+\}',body[i:])
            if nested:env+=1 if nested[1]=='begin' else -1;i+=nested.end();continue
            i+=2;continue
        depth+=int(body[i]=='{')-int(body[i]=='}');i+=1
    out.append(body[start:]);return out


class Preprocessor:
    def __init__(self,chapter,number,registry):
        self.chapter=chapter;self.number=number;self.equation=0;self.theorem=0;self.registry=registry
        self.section=[0,0,0];self.figure_ids=[]
    def file(self,path):
        text=active(path.read_text())
        # Replace figures before expanding includes. Their contents remain available
        # in the approved asset, and never go through a lossy HTML TeX conversion.
        fs=[f for f in figures(path,self.chapter) if f['source']==str(path.relative_to(ROOT))]
        for f in reversed(fs):
            key=next((x['id'] for x in self.registry if x.get('book_label') and x['book_label']==f['book_label']),None)
            if key is None:
                paths=[g['reference'] for g in f['assets']]
                key=next((x['id'] for x in self.registry if x['chapter']==self.chapter and any(x['folder'] in p for p in paths)),None)
            if key is None:raise RuntimeError(f'Unmapped figure: {f["source"]}:{f["line"]}')
            self.figure_ids.append(key)
            text=text[:f['start']]+'\n\\webfigure{'+key+'}\n'+text[f['end']:]
        def include(name):
            p=ROOT/name
            if not p.suffix:p=p.with_suffix('.tex')
            return self.file(p)
        text=replace_commands(text,'input',include)
        for g in reversed(graphics(text)):
            item=next((x for x in self.registry if x['chapter']==self.chapter and x['folder'] in g['reference']),None)
            if item is None:raise RuntimeError('Unmapped inline illustration '+g['reference'])
            self.figure_ids.append(item['id'])
            text=text[:g['start']]+'\n\\webfigure{'+item['id']+'}\n'+text[g['end']:]
        return text
    def prepare(self,text):
        # Tables and lists retain their contents; page-layout directives do not
        # constrain a reflowable edition.
        text=re.sub(r'\\(?:begin|end)\{samepage\}','',text)
        text=text.replace(r'\begin{rs}',r'\begin{itemize}').replace(r'\end{rs}',r'\end{itemize}')
        text=replace_commands(text,'chapteroverview',lambda x:'\n\\weboverview\n\n'+x+'\n\n\\webendoverview\n')
        text=replace_commands(text,'texorpdfstring',lambda x,y:x,2)
        # Emit exact counters before Pandoc discards equation/align distinctions.
        pattern=re.compile(r'\\(eqml|eqm|eql|eq)\b|\\begin\{(equation\*?|align\*?)\}')
        result=[];pos=0
        for m in pattern.finditer(text):
            if m.start()<pos:continue
            if m[1]:body,end=braced(text,m.end());kind={'eq':'equation*','eql':'equation','eqm':'align*','eqml':'align'}[m[1]]
            else:
                kind=m[2];closing='\\end{'+kind+'}';end=text.index(closing,m.end());body=text[m.end():end];end+=len(closing)
            isalign=kind.startswith('align');numbered=not kind.endswith('*');parts=split_rows(body) if isalign and numbered else [body]
            rendered=[]
            for part in parts:
                if not part.strip():continue
                tagged=''
                manual=re.search(r'\\tag\*?\s*\{',part)
                if manual:
                    tag,tagend=braced(part,manual.end()-1)
                    tagged=r'\tag{'+tag+'}'
                    part=part[:manual.start()]+part[tagend:]
                elif numbered and not re.search(r'\\(?:nonumber|notag)\b',part):
                    self.equation+=1;number=f'{self.number}.{self.equation}';tagged=r'\tag{'+number+'}'
                    for lab in re.findall(r'\\label\{([^}]+)\}',part):
                        if lab in LABELS and LABELS[lab][0]!=number:raise RuntimeError(f'Equation counter drift: {lab}: {number} vs {LABELS[lab][0]}')
                    REPORT['equations'].append({'chapter':self.chapter,'number':number,'labels':re.findall(r'\\label\{([^}]+)\}',part)})
                part=re.sub(r'\\(?:nonumber|notag)\b','',part)
                if isalign:part='\\begin{aligned}\n'+part+'\n\\end{aligned}'
                rendered.append('\\[\n'+part+'\n'+tagged+'\n\\]')
            result.extend([text[pos:m.start()],'\n\n'.join(rendered)]);pos=end
        result.append(text[pos:]);text=''.join(result)
        # Pandoc omits optional titles on custom theorem environments. Preserve
        # them explicitly, and keep the shared theorem counter used by the PDF.
        pattern=re.compile(r'\\(begin|end)\{('+ '|'.join([*THEOREMS,'proof'])+r')\}(?:\[([^\]]*)\])?')
        def theorem(m):
            if m[1]=='end':return '\n\n\\webendbox\n\n'
            typ=m[2];title=m[3] or '';number=''
            if typ in THEOREMS and typ!='ass':self.theorem+=1;number=f'{self.number}.{self.theorem}'
            name='Proof' if typ=='proof' else THEOREMS[typ]+(' '+number if number else '')
            if title:name=(title if typ=='proof' else name+' ('+title+')')
            return '\n\n\\webbeginbox{'+typ+'}{'+name+'}\n\n'
        text=pattern.sub(theorem,text)
        return text


def tex_math(text):
    text=re.sub(r'\\label\{[^}]*\}', '', text)
    text=text.replace(r'\mathds',r'\mathbb').replace(r'\llbracket',r'\lbrack\!\lbrack').replace(r'\rrbracket',r'\rbrack\!\rbrack')
    text=text.replace(r'\upshape','').replace(r'\tiny','').replace(r'\mbox',r'\text')
    # A TeX accent consumes a macro token before expansion. Pandoc expands
    # \hat\RR to \hat\mathbb{R}; group the expanded argument for browser TeX.
    pattern=re.compile(r'\\(hat|bar|tilde|vec|dot|ddot|widehat|widetilde|overline|underline)\s*(\\(?:mathbb|mathcal|mathscr|mathrm|mathbf|mathsf)\s*\{)')
    pieces=[];pos=0
    for m in pattern.finditer(text):
        if m.start()<pos:continue
        _,end=braced(text,m.end()-1)
        pieces.extend([text[pos:m.start()],'\\'+m[1]+'{'+text[m.start(2):end]+'}']);pos=end
    pieces.append(text[pos:]);text=''.join(pieces)
    return text.strip()


class Renderer:
    def __init__(self,registry):
        self.registry={x['id']:x for x in registry};self.section=[0,0,0];self.boxdepth=0;self.seen_labels=set();self.nodecounts=collections.Counter()
    def anchor(self,label):
        label=web_id(ALIASES.get(label,label))
        if label in self.seen_labels:return ''
        self.seen_labels.add(label);return '\n\n('+label+')=\n'
    def raw(self,s):
        ref=re.match(r'\\(ref|eqref|pageref|label)\{([^}]+)\}',s)
        if ref:
            return self.anchor(ref[2]) if ref[1]=='label' else self.linkref(ref[2],ref[1]=='eqref',ref[1]=='pageref')
        if s.strip() in [r'\textquotesingle',r'\textregistered',r'\quad']:
            return {r'\textquotesingle':"'",r'\textregistered':'®',r'\quad':' '}[s.strip()]
        if s.startswith(r'\webfigure'):
            key,_=braced(s,len(r'\webfigure'));f=self.registry[key];num=f.get('figure_number');label=f.get('book_label') or 'fig-'+key
            REPORT['figures'].append(key)
            cap=self.inline_blocks(pandoc(f['caption_source']))
            anchor=self.anchor(label)
            prefix='Figure '+num+'. ' if num else 'Illustration. '
            return anchor+'\n:::{div}\n:class: book-figure\n\n'+f'<iframe class="book-live-frame" title="{html.escape(prefix+plain(pandoc(f["caption_source"]))[:140],quote=True)}" src="../live/panels/{key}.html" loading="lazy" style="width:100%;height:560px;border:0;display:block;"></iframe>\n\n'+f'**{prefix.strip()}** {cap}\n\n:::\n\n'
        if s.startswith(r'\webbeginbox'):
            typ,p=braced(s,len(r'\webbeginbox'));name,_=braced(s,p);self.boxdepth+=1
            # Four-colon fences allow figures and ordinary directives inside proofs.
            fence=':'*(4+self.boxdepth)
            return '\n\n'+fence+'{admonition} '+self.inline_blocks(pandoc(name))+'\n:class: '+('proof' if typ=='proof' else 'theorem')+'\n\n'
        if s.startswith(r'\webendbox'):
            fence=':'*(4+self.boxdepth);self.boxdepth-=1;return '\n\n'+fence+'\n\n'
        if s==r'\weboverview':return '\n\n:::{div}\n:class: chapter-overview\n\n'
        if s==r'\webendoverview':return '\n\n:::\n\n'
        if re.match(r'\\(?:medskip|smallskip|bigskip|noindent|centering|newpage|clearpage|vfill|hfill|par|relax|protect|BookContinuedFigure)\b',s):return '\n'
        if re.match(r'\\(?:newcommand|renewcommand|DeclareMathOperator|def)\b',s):return ''
        if s in [r'\begin{center}',r'\end{center}',r'\begin{abstract}',r'\end{abstract}']:return '\n'
        if re.match(r'\\(?:vspace|hspace|setlength|enlargethispage|bibliographystyle)\*?(?:\[.*?\])?\{',s):return ''
        REPORT['unknown_tex'].append({'chapter':PAGE,'tex':s});return s
    def linkref(self,key,eq=False,page=False):
        REPORT['references'].append({'chapter':PAGE,'target':key})
        if key not in LABELS:raise RuntimeError('Missing publication reference '+key)
        number=LABELS[key][1 if page else 0];canonical=ALIASES.get(key,key);owner=OWNERS[canonical];target=(owner+'.md#'+web_id(canonical)) if owner!=PAGE else '#'+web_id(canonical)
        if number.startswith('{') and number.endswith('}'):number=number[1:-1]
        text=('('+number+')') if eq else number
        return '['+text+']('+target+')'
    def inline(self,nodes):
        out=[]
        for n in nodes:
            t=n['t'];c=n.get('c');self.nodecounts[t]+=1
            if t=='Str':out.append(c.replace('[',r'\[').replace(']',r'\]').replace('*',r'\*'))
            elif t in ('Space','SoftBreak'):out.append(' ')
            elif t=='LineBreak':out.append('  \n')
            elif t in ('Emph','Strong','SmallCaps','Strikeout'):
                mark={'Emph':'*','Strong':'**','SmallCaps':'','Strikeout':'~~'}[t];out.append(mark+self.inline(c)+mark)
            elif t in ('Superscript','Subscript'):out.append('$'+('^' if t=='Superscript' else '_')+'{'+self.inline(c)+'}$')
            elif t=='Quoted':out.append('“'+self.inline(c[1])+'”')
            elif t=='Math':
                labels=re.findall(r'\\label\{([^}]+)\}',c[1]);math=tex_math(c[1])
                if c[0]['t']=='DisplayMath':out.append(''.join(self.anchor(l) for l in labels)+'\n\n```{math}\n'+math+'\n```\n\n')
                else:out.append('$'+math+'$')
            elif t=='RawInline':out.append(self.raw(c[1]))
            elif t=='Code':out.append('`'+c[1]+'`')
            elif t=='Span':
                if c[0][0]:out.append(self.anchor(c[0][0]))
                out.append(self.inline(c[1]))
            elif t=='Link':
                attr=dict(c[0][2]);typ=attr.get('reference-type')
                if typ:out.append(self.linkref(attr['reference'],typ=='eqref',typ=='pageref'))
                else:out.append('['+self.inline(c[1])+']('+c[2][0]+')')
            elif t=='Cite':
                refs=[]
                for cite in c[0]:
                    key=cite['citationId'];num=BIB.get(key)
                    if num is None:raise RuntimeError('Unresolved citation '+key)
                    refs.append('['+num+'](../references.md#bib-'+key+')')
                out.append('['+', '.join(refs)+']')
            elif t=='Note':out.append(' ('+self.inline_blocks(c).strip()+')')
            else:raise RuntimeError(f'Unimplemented inline {t}: {str(c)[:250]}')
        return ''.join(out)
    def inline_blocks(self,nodes):return self.blocks(nodes).strip()
    def blocks(self,nodes):
        out=[]
        for n in nodes:
            t=n['t'];c=n.get('c');self.nodecounts[t]+=1
            if t in ('Para','Plain'):out.append(self.inline(c)+'\n\n')
            elif t=='Header':
                level,attr,content=c;title=self.inline(content)
                if level==1:continue
                if level>=4:
                    if attr[0] in LABELS:out.append(self.anchor(attr[0]))
                    out.append('**'+title+'**\n\n');continue
                if level in (2,3):
                    self.section[level-2]+=1
                    if level==2:self.section[1]=0
                    prefix=str(CHNO)+'.'+'.'.join(str(k) for k in self.section[:level-1]);title=prefix+' '+title
                if attr[0]:out.append(self.anchor(attr[0] if attr[0] in LABELS else PAGE+'-'+attr[0]))
                out.append('#'*min(level,5)+' '+title+'\n\n')
            elif t=='RawBlock':out.append(self.raw(c[1])+'\n\n')
            elif t=='Div':
                if c[0][0]:out.append(self.anchor(c[0][0]))
                out.append(self.blocks(c[1]))
            elif t=='BulletList':
                for item in c:out.append('- '+self.blocks(item).strip().replace('\n','\n  ')+'\n')
                out.append('\n')
            elif t=='OrderedList':
                for i,item in enumerate(c[1],c[0][0]):out.append(str(i)+'. '+self.blocks(item).strip().replace('\n','\n   ')+'\n')
                out.append('\n')
            elif t=='DefinitionList':
                for term,defs in c:out.append('**'+self.inline(term)+'**\n\n'+''.join(self.blocks(d) for d in defs))
            elif t=='BlockQuote':out.append('> '+self.blocks(c).strip().replace('\n','\n> ')+'\n\n')
            elif t=='CodeBlock':out.append('```'+(' '.join(c[0][1]))+'\n'+c[1]+'\n```\n\n')
            elif t=='HorizontalRule':out.append('---\n\n')
            elif t=='Table':
                # Delegate the rare text table to Pandoc's HTML writer; it retains
                # row and column spans rather than flattening them into prose.
                obj={'pandoc-api-version':[1,23,1],'meta':{},'blocks':[n]}
                val=subprocess.run(['pandoc','-f','json','-t','html','--mathjax'],input=json.dumps(obj),text=True,capture_output=True,check=True).stdout
                out.append(val+'\n\n')
            elif t=='Null':pass
            else:raise RuntimeError(f'Unimplemented block {t}: {str(c)[:250]}')
        return re.sub(r'\n{4,}','\n\n\n',''.join(out))


def plain(nodes):
    def p(v):
        if isinstance(v,list):return ''.join(p(x) for x in v)
        if isinstance(v,dict):
            if 't' not in v:return ''
            if v['t']=='Str':return v['c']
            if v['t'] in ['Space','SoftBreak']:return ' '
            if v['t']=='Math':return v['c'][1]
            return p(v.get('c',[]))
        return ''
    return p(nodes)


def registry():
    old=json.loads((ROOT/'docs/figures/figure-project.json').read_text())
    live=[x for x in old if x['chapter']!='cover' and x['review_state'] not in ['removed','merged']]
    # Use CURRENT captions and numbers, not numbers from the previous review.
    current=[]
    for ch in chapter_list():current.extend(figures(ROOT/'chapters-tex'/f'{ch}.tex',ch))
    used=set();out=[]
    for f in current:
        snippet=active((ROOT/f['source']).read_text())[f['start']:f['end']]
        for label in re.findall(r'\\label\{([^}]+)\}',snippet):
            if label!=f['book_label'] and f['book_label']:ALIASES[label]=f['book_label']
        item=next((x for x in live if x.get('book_label') and x['book_label']==f['book_label'] and x['id'] not in used),None)
        if item is None:item=next((x for x in live if x['id'] not in used and x['chapter']==f['chapter'] and any(x['folder'] in a['reference'] for a in f['assets'])),None)
        if item is None:raise RuntimeError(f'No manifest figure {f["source"]}:{f["line"]}')
        item=dict(item);used.add(item['id']);item['caption_source']=f['caption_source'] or item['caption_source'];item['book_label']=f['book_label']
        if f['book_label'] in LABELS:item['figure_number']=LABELS[f['book_label']][0];item['book_page']=LABELS[f['book_label']][1]
        item['current_assets']=[a['path'] for a in f['assets']]
        out.append(item)
    for item in live:
        if item['id'] in used:continue
        item=dict(item)
        raw=active((ROOT/'chapters-tex'/f'{item["chapter"]}.tex').read_text())
        item['current_assets']=[g['path'] for g in graphics(raw) if item['folder'] in g['reference']]
        if not item['current_assets']:raise RuntimeError('Missing unnumbered illustration '+item['id'])
        out.append(item)
    return out


def main():
    global LABELS,OWNERS,BIB,PAGE,CHNO
    data=references();LABELS=data['labels'];OWNERS=data['owners'];BIB=data['bibliography']
    by_anchor={}
    for key,fields in LABELS.items():
        if len(fields)>3 and fields[3]:
            anchor=(OWNERS[key],fields[3])
            if anchor in by_anchor:ALIASES[key]=by_anchor[anchor]
            else:by_anchor[anchor]=key
    reg=registry();(WEB/'chapters').mkdir(parents=True,exist_ok=True);(WEB/'live').mkdir(exist_ok=True)
    chapters=[];coverage=json.loads((WEB/'interactions.json').read_text()) if (WEB/'interactions.json').exists() else {}
    public=[]
    for f in reg:
        obj={k:f.get(k) for k in ['id','chapter','chapter_title','figure_number','book_page','book_label','caption_source','folder','engine','review_state','current_assets','asset_directory','code_directory']}
        obj['interaction']=coverage.get(f['id']);public.append(obj)
    (WEB/'live/figures.json').write_text(json.dumps(public,ensure_ascii=False,indent=2)+'\n')
    for CHNO,PAGE in enumerate(chapter_list(),1):
        pre=Preprocessor(PAGE,CHNO,reg);path=ROOT/'chapters-tex'/f'{PAGE}.tex';source=pre.prepare(pre.file(path));ast=pandoc(source)
        renderer=Renderer(reg);header=next(n for n in ast if n['t']=='Header' and n['c'][0]==1);title=renderer.inline(header['c'][2]);anchor=header['c'][1][0]
        intro='---\n'+f'title: {json.dumps(str(CHNO)+". "+title)}\n'+'numbering:\n  headings: false\n  math: false\n  figures: false\n---\n\n'
        intro+=renderer.anchor(anchor) if anchor else ''
        intro+=f'[Chapter PDF](../assets/pdfs/chapters-pdf/{PAGE}.pdf) · [Complete book](../assets/pdfs/FundationsDataScience.pdf)\n\n'
        body=renderer.blocks(ast)
        if renderer.boxdepth:raise RuntimeError('Unclosed theorem/proof in '+PAGE)
        (WEB/'chapters'/f'{PAGE}.md').write_text(intro+body)
        chapters.append({'file':f'chapters/{PAGE}.md'})
        REPORT['chapters'].append({'chapter':PAGE,'title':title,'number':CHNO,'source_sha256':hashlib.sha256(source.encode()).hexdigest(),'equations':pre.equation,'theorems':pre.theorem,'figures':len(pre.figure_ids),'words':len(plain(ast).split()),'nodes':dict(renderer.nodecounts)})
    PAGE='references';CHNO=0
    bbl=ROOT/'build/book/FundationsDataScience.bbl';saved=WEB/'references-source.tex'
    if bbl.exists():shutil.copyfile(bbl,saved)
    text=saved.read_text();items=list(re.finditer(r'\\bibitem(?:\[[^\]]*\])?\{([^}]+)\}',text));render=Renderer(reg);parts=['# References\n\n']
    for i,m in enumerate(items):
        body=text[m.end():items[i+1].start() if i+1<len(items) else text.index(r'\end{thebibliography}')]
        body=body.replace(r'\newblock',' ')
        parts.append(f'(bib-{m[1]})=\n\n**[{BIB[m[1]]}]** '+render.blocks(pandoc(body)))
    (WEB/'references.md').write_text(''.join(parts))
    cover='# Mathematical Foundations of Data Sciences\n\n'
    cover+='Sampling, representation, information, and learning: read the mathematics, then explore how its conclusions change as you vary the parameters. Each experiment sits beside its corresponding book figure and retains the publication’s numbering.\n\n'
    cover+='[Download the complete book](assets/pdfs/FundationsDataScience.pdf)\n\n'
    cover+=':::{figure} assets/cover.jpg\n:alt: Original wave illustration from the book cover\n:width: 85%\n:::\n\n'
    preface=active((ROOT/'chapters-tex/abstract.tex').read_text());preface=re.sub(r'\\chapter\*\{[^}]*\}|\\addcontentsline\{[^}]*\}\{[^}]*\}\{[^}]*\}|\\tableofcontents','',preface)
    cover+='## Preface\n\n'+Renderer(reg).blocks(pandoc(preface))
    cover+='## Chapters\n\n'+''.join(f'{i["number"]}. [{i["title"]}](chapters/{i["chapter"]}.md)\n' for i in REPORT['chapters'])
    (WEB/'index.md').write_text(cover)
    cfg={'version':1,'project':{'id':'mathematical-foundations-of-data-sciences','title':'Mathematical Foundations of Data Sciences','authors':[{'name':'Gabriel Peyré','affiliations':[{'name':'CNRS & DMA, École Normale Supérieure'}]}],'description':'The complete book with interactive experiments in sampling, wavelets, information, inverse problems, and learning.','toc':[{'file':'index.md'},*chapters,{'file':'references.md'}],'static_files':['live','assets','book-frame.js'],'numbering':{'headings':False,'math':False,'figures':False}},'site':{'template':'book-theme','options':{'logo_text':'Foundations of Data Sciences','style':'book-theme.css','hide_footer_links':False,'outline_maxdepth':3},'parts':{'primary_sidebar_footer':''}}}
    (WEB/'myst.yml').write_text(json.dumps(cfg,ensure_ascii=False,indent=2)+'\n')
    (WEB/'conversion-report.json').write_text(json.dumps(REPORT,indent=2,ensure_ascii=False)+'\n')
    print(f'Converted {len(chapters)} chapters, {len(REPORT["equations"])} numbered equations, {len(REPORT["figures"])} figure displays.')
    if REPORT['unknown_tex']:
        for item in REPORT['unknown_tex'][:30]:print('Unhandled:',item)
        raise SystemExit(f'{len(REPORT["unknown_tex"])} unhandled TeX constructs; inspect conversion-report.json')

if __name__=='__main__':main()
