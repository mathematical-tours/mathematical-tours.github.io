# Website refresh — September 2026

The website now presents the three books as its main collection, followed by Numerical Tours, Mathematical Nexus, the algorithms directory, course slides, and the ten-lecture machine-learning course.

## Main changes

- Shared navigation, footer, metadata, typography, and colorful artwork banners across all seven main pages. Book pages keep their full-book downloads, chapter directories, and edition-specific resources.
- Responsive book cards explain each book’s subject, intended background, and available formats. The introduction remains available in English and French; DAFT’s relationship to the original French Ellipses book remains explicit.
- Renamed the course to **Mathematics of Machine Learning**, preserving `/maths-ia-course/`. Its ten lectures have a contents menu, topic lists, practical materials, and further reading. The page distinguishes the nine available transcripts from the tenth lecture’s OT4ML materials.
- Rebuilt Algorithms as a searchable, alphabetized directory of 216 entries in ten subjects. Replaced the previous chronology and single-person attributions with linked references, avoiding misleading claims about dates or sole inventors. Corrected names and conflations including Lempel–Ziv versus LZMA, bagging versus bootstrap, Burer–Monteiro, and the distinction between GMRES and BiCGSTAB.
- Reorganized all seventeen slide decks into two ordered collections with short descriptions.
- Repaired moved optimization-note paths, retired notebook hosts, the automatic-differentiation slide link, and several research references. The old Tweets pages forward to Mathematical Nexus.
- Replaced the old Bootstrap/jQuery navigation and obsolete analytics scripts with lightweight shared HTML, CSS, and progressive enhancement. Navigation and catalogues remain usable without JavaScript. Added a consistent 404 page and an SVG favicon.
- Excluded development and legacy source folders from the website build while retaining the public books, chapters, course materials, Python archive, and DAFT figure comparator. In particular, the MFDS figure-development archive no longer adds 231 MB to the published site.
- Updated the introduction’s website templates so regenerating its page preserves the shared design.

## Verification

- Complete Jekyll build with safe mode succeeds.
- `tools/check-site.py` checks eleven public routes and 95 distinct local targets: no missing pages, assets, chapter PDFs, or static anchors; no duplicate IDs or missing image alternatives.
- Seven main pages reviewed at widths of 1360, 768, 390, and 320 pixels, with no horizontal overflow or browser script errors. Search, subject selection, empty results, reset controls, mobile navigation, Escape behavior, and the skip link were exercised. The chapter and algorithm directories were also checked with JavaScript disabled.
- Full-book and standalone-chapter PDF metadata checked. All eight DAFT links open the intended chapter starts; all nineteen Foundations chapter lengths match their PDFs.
- The interactive edition’s home and Shannon chapter load from the complete publishing build without failed local requests.
- External course, slide, and algorithm links audited. Broken references were replaced; some publisher sites, IEEE, and the OpenAI tokenizer limit automated requests, so HTTP access alone is not treated as proof that a reference is unavailable.
- The existing book archives were verified before publication: 856 DAFT files and 335 introductory-book files match their original SHA-256 manifests.

The book revisions have their own correction logs and build/test reports inside their source folders. Website source locations and local build instructions are in the root README.
