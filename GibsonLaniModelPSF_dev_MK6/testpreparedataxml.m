clc
clear

docNode = com.mathworks.xml.XMLUtils.createDocument('VutaraDataConfiguration');
%%
VutaraDataConfiguration = docNode.getDocumentElement;
VutaraDataConfiguration.setAttribute('version','1.0');
product=docNode.createElement('VutaraDataConfigurationitem');
product.setAttribute('target','avs');
VutaraDataConfiguration.appendChild(product);

xmlwrite('infoUAT.xml',docNode);
type('infoUAT.xml');
%%
docNode = com.mathworks.xml.XMLUtils.createDocument('VDC');
VDC = docNode.getDocumentElement;
VDC.setAttribute('version','1.0')

