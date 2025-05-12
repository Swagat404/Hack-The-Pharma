# 🚀 Hack-The-Pharma!

Our team, **Hack-The-Pharma**, developed a data-driven solution to **Challenge 1: Predict Kinase Selectivity Using Machine Learning for PharmaHacks 2025**. [Challenge Repository](https://github.com/devpatel30/molecular_forecaster). Our project aims to address critical challenges in pharmaceutical data analysis and molecular property prediction, using Python and R. In drug development 💊, predicting molecular properties is crucial for screening potential drug candidates efficiently. Our solution leverages machine learning techniques to forecast kinase-inhibitor interactions from molecular properties and chemical structures.

---

## 📦 Installation and Setup

1. **Clone the Repository**

```bash
git clone https://github.com/Swagat404/Hack-The-Pharma.git
cd Hack-The-Pharma
```

2. **Create a Virtual Environment (Optional)**

```bash
python -m venv venv
source venv/bin/activate   # On macOS/Linux
.\venv\Scripts\activate    # On Windows
```

3. **Install Python Dependencies**

```bash
pip install -r requirements.txt
```

4. **Install R Packages**
   From your R console:

```r
install.packages(c("dplyr", "ggplot2", "caret", "randomForest"))
```

---

## ▶️ Running the Project

To run the main pipeline:

```bash
python main.py
```

To execute R scripts:

```bash
Rscript script_name.R
```

## 🚧 Future Improvements

* Enhanced feature selection for better model accuracy
* Integration with external chemical databases to enrich available input data
* Better script organization for testing calls
