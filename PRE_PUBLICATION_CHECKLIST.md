# MetaboGraph Pre-Publication Checklist

## 📋 Essential Items

### ✅ Code & Documentation (COMPLETED)
- [x] All code files present and functional
- [x] User guides created for all tabs
- [x] README.md comprehensive and accurate
- [x] CONTRIBUTING.md accurate
- [x] LICENSE file with database licenses
- [x] Python version requirements correct (3.8+)

---

## 🔴 REQUIRED Before Publishing

### 1. Personal Information & Branding
- [ ] **Update README.md placeholders:**
  - [ ] Replace `[Your Name]` with actual name
  - [ ] Replace `your.email@example.com` with real email
  - [ ] Replace `yourusername` in GitHub URLs
  - [ ] Replace `support@metabograph.org` with real contact
  - [ ] Update `@software{metabograph2024, author = {Your Name}` in citation

- [ ] **Update setup.py:**
  - [ ] Line 18: `author="Your Name"` → Real name
  - [ ] Line 19: `author_email="your.email@example.com"` → Real email
  - [ ] Line 23: GitHub URL

- [ ] **Update LICENSE:**
  - [ ] Line 3: `Copyright (c) 2024 [Your Name/Organization]` → Real info

### 2. Screenshots & Images
**Location:** `docs/images/`

**Required Screenshots:**
- [ ] **Main window** - All tabs visible
- [ ] **Data Cleaning tab** - Show interface with sample data loaded
- [ ] **ID Annotation tab** - Show annotation results
- [ ] **Pathway Annotation tab** - Show pathway mapping progress
- [ ] **Pathway Network tab** - Show enrichment results and network
- [ ] **Comparative Analysis tab** - Show Z-score heatmap
- [ ] **Multi-Omics Integration tab** - Show merge interface
- [ ] **Database Setup tab** - Show database status

**Format:** PNG or JPG, 1200-1600px width recommended

**Add to README.md:** Replace line 65 `*[Add screenshots of your application here]*`

### 3. Example/Test Data Files
**Location:** Create `examples/` or `test_data/` folder

**Recommended Files from Your Manuscript:**
- [ ] **Sample cleaned data** (Excel)
  - Small dataset (50-100 metabolites)
  - Both positive and negative mode sheets
  - Ready for ID annotation
  - File: `example_cleaned_metabolites.xlsx`

- [ ] **Sample ID-annotated data** (Excel)
  - With HMDB, KEGG, PubChem IDs
  - Ready for pathway annotation
  - File: `example_id_annotated.xlsx`

- [ ] **Sample pathway-annotated data** (Excel)
  - With pathway columns
  - Ready for network analysis
  - File: `example_pathway_annotated.xlsx`

- [ ] **Sample comparative data** (2-3 Excel files)
  - For testing comparative analysis
  - Each with "Selected_Pathways" sheet
  - Files: `dataset1_pathways.xlsx`, `dataset2_pathways.xlsx`

- [ ] **README in examples folder**
  - File: `examples/README.md`
  - Describe each file and how to use it

### 4. Test Documentation (RECOMMENDED)
**Location:** Create `TESTING.md` or `docs/testing.md`

**Should Include:**
- [ ] **Manual Testing Results**
  - Each tab tested with example data
  - Expected vs. actual results
  - Known issues/limitations

- [ ] **Test Cases:**
  ```markdown
  ## Data Cleaning Tab
  - ✅ Load Compound Discoverer export
  - ✅ Process positive mode
  - ✅ Process negative mode
  - ✅ Reference ion filtering works
  - ✅ Export cleaned data
  
  ## ID Annotation Tab
  - ✅ Annotate metabolites (HMDB)
  - ✅ API queries work (PubChem)
  - ✅ Lipid annotation works
  - ✅ Export annotated data
  
  [Continue for each tab...]
  ```

### 5. Database Files Instructions
- [ ] **Create database download guide:**
  - File: `docs/DATABASE_SETUP_GUIDE.md`
  - Step-by-step with screenshots
  - Links to all database sources
  - Processing instructions

- [ ] **Note in README:** Users must download databases separately (due to size)

### 6. Distribution Files
**For Zenodo/OneDrive Upload:**

- [ ] **Create release package:**
  ```
  MetaboGraph_v1.0.1.zip containing:
  ├── MetaboliteAnnotationTool.exe
  ├── _internal/ (PyInstaller files)
  ├── README.txt (quick start)
  ├── LICENSE.txt
  └── Databases/ (empty folder with README)
  ```

- [ ] **Create README.txt for release:**
  ```
  MetaboGraph v1.0.1
  ==================
  
  QUICK START:
  1. Extract all files
  2. Run MetaboliteAnnotationTool.exe
  3. Go to Database Setup tab
  4. Download and process databases
  5. Start analysis!
  
  Full Documentation: See online documentation
  Support: your.email@example.com
  ```

### 7. Additional Documentation Files
- [ ] **CONTRIBUTORS.md** - List contributors (if any)
- [ ] **CHANGELOG.md** - Version history
- [ ] **CITATION.cff** - Structured citation file (for GitHub)
- [ ] **CODE_OF_CONDUCT.md** - Community guidelines

### 8. Final Code Checks
- [ ] **Remove debug code:**
  - [ ] Check for `print()` statements
  - [ ] Remove test/development comments
  - [ ] Remove unused imports

- [ ] **Verify all paths:**
  - [ ] Database paths work correctly
  - [ ] Output paths are configurable
  - [ ] No hardcoded absolute paths

- [ ] **Test executable:**
  - [ ] Test on clean Windows machine (without Python)
  - [ ] Verify all tabs open
  - [ ] Test basic workflow end-to-end

### 9. Legal & Ethics
- [ ] **Verify database licenses**
  - Confirm you can redistribute (most are OK)
  - Cite all databases properly

- [ ] **Check dependencies**
  - All packages in requirements.txt are open source
  - No GPL conflicts with MIT license

---

## 🟡 RECOMMENDED (Enhance Quality)

### Additional Features
- [ ] **Video tutorial** (YouTube, 5-10 min walkthrough)
- [ ] **FAQ.md** - Common questions
- [ ] **Troubleshooting.md** - Common issues
- [ ] **Performance benchmarks** - Processing times
- [ ] **Comparison with other tools** (MetaboAnalyst, etc.)

### Code Quality
- [ ] Add docstrings to all functions
- [ ] Add type hints
- [ ] Create unit tests (optional but professional)
- [ ] Run linter (pylint, flake8)

### Community
- [ ] Create GitHub Discussions/Issues templates
- [ ] Add badges to README (downloads, version, etc.)
- [ ] Create social media presence (Twitter, ResearchGate)

---

## 📦 Publication Platforms Checklist

### Zenodo Upload
- [ ] Create account
- [ ] Upload release package
- [ ] Add metadata (title, authors, description, keywords)
- [ ] Select license (MIT)
- [ ] Add related publications (manuscript DOI)
- [ ] Publish and get DOI
- [ ] Update README with Zenodo DOI link

### Alternative: OneDrive
- [ ] Upload release package
- [ ] Create shareable link
- [ ] Test link in incognito browser
- [ ] Update README with download link

---

## ✅ Final Verification

### Before Going Live:
1. [ ] Download your own release package
2. [ ] Extract on fresh computer (or VM)
3. [ ] Follow your own README instructions
4. [ ] Complete one full analysis workflow
5. [ ] Verify all documentation links work
6. [ ] Check all external URLs are valid

### Post-Publication:
- [ ] Share on social media/academic networks
- [ ] Email collaborators
- [ ] Add to CV/portfolio
- [ ] Submit manuscript citing software

---

## 📊 Priority Ranking

**CRITICAL (Must Do):**
1. Update all personal information (names, emails, URLs)
2. Add screenshots to README
3. Create example data files
4. Test executable on clean machine
5. Upload to distribution platform

**HIGH PRIORITY (Should Do):**
6. Create TESTING.md with test results
7. Add examples/README.md
8. Create release package with quick start
9. Verify all documentation

**NICE TO HAVE (Optional):**
10. Video tutorial
11. Unit tests
12. Additional documentation files

---

## 🎯 Time Estimates

- **Critical items:** 4-6 hours
- **High priority:** 3-4 hours
- **Nice to have:** 8-12 hours

**Total for publication-ready:** ~7-10 hours

---

## 📝 Notes

- Your manuscript data files are PERFECT for examples
- Include small subset (not full dataset) for testing
- Add DOI to manuscript in README once published
- Consider embargo timing if coordinating with paper publication
