% Build the complete book and all active chapters with the shared build tool.
bookRoot = fileparts(fileparts(mfilename('fullpath')));
previousDirectory = pwd;
restoreDirectory = onCleanup(@() cd(previousDirectory));
cd(bookRoot);
status = system('python3 scripts/build_book.py');
assert(status == 0, 'Book compilation failed; inspect build/ for diagnostics.');
